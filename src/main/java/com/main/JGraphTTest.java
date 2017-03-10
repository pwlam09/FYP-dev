package com.main;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.geom.Rectangle2D;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Properties;
import java.util.Set;

import javax.swing.JApplet;
import javax.swing.JFrame;

import org.jgraph.JGraph;
import org.jgraph.graph.AttributeMap;
import org.jgraph.graph.DefaultGraphCell;
import org.jgraph.graph.GraphConstants;
import org.jgrapht.DirectedGraph;
import org.jgrapht.ListenableGraph;
import org.jgrapht.WeightedGraph;
import org.jgrapht.ext.JGraphModelAdapter;
import org.jgrapht.graph.DefaultDirectedWeightedGraph;
import org.jgrapht.graph.DefaultListenableGraph;
import org.jgrapht.graph.DefaultWeightedEdge;

import edu.stanford.nlp.ling.CoreAnnotations.PartOfSpeechAnnotation;
import edu.stanford.nlp.ling.CoreAnnotations.SentencesAnnotation;
import edu.stanford.nlp.ling.CoreAnnotations.TextAnnotation;
import edu.stanford.nlp.ling.CoreAnnotations.TokensAnnotation;
import edu.stanford.nlp.ling.CoreLabel;
import edu.stanford.nlp.pipeline.Annotation;
import edu.stanford.nlp.pipeline.StanfordCoreNLP;
import edu.stanford.nlp.util.CoreMap;
import edu.stanford.nlp.util.StringUtils;

public class JGraphTTest extends JApplet {
	private static JGraphTTest instance = new JGraphTTest();
    private static final long serialVersionUID = 3256444702936019250L;
    private static final Color DEFAULT_BG_COLOR = Color.decode("#FAFBFF");
    private static final Dimension DEFAULT_SIZE = new Dimension(530, 320);
    
    // directed graph
	private static ListenableDirectedWeightedGraph<WordPosTuple, CustomDirectedWeightedEdge> wordGraph = 
			new ListenableDirectedWeightedGraph<>(CustomDirectedWeightedEdge.class);
	private static WordPosTuple startWord;
	private static WordPosTuple endWord;
	
	// stopwords will be loaded and concatenated with the separator
	private static String stopWordSeparator = ",";
	private static String linkedStopwords;
	
	private static String punct_tag = "PUNCT";
	
	// minimum number of words for compression
	private static int nb_words = 8;
	
	// at least one verb for each compression
	private static String[] verbs = {"VB", "VBD", "VBP", "VBZ", "VH", "VHD", "VHP", "VV", "VVD", "VVP", "VVZ"};

	
	private static HashMap<WordPosTuple, Integer> term_freq = new HashMap<>();
	
    private JGraphModelAdapter<WordPosTuple, CustomDirectedWeightedEdge> jgAdapter;

    private JGraphTTest() {
    	
    }
    
    /**
     * An alternative starting point for this demo, to also allow running this applet as an
     * application.
     *
     * @param args ignored.
     */
    public static void main(String[] args)
    {
    	JGraphTTest applet = new JGraphTTest();
    	loadStopwordsList();
        applet.init();

        JFrame frame = new JFrame();
        frame.getContentPane().add(applet);
        frame.setTitle("JGraphT Adapter to JGraph Demo");
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.pack();
        frame.setVisible(true);
    }

    private static void loadStopwordsList() {
    	if (linkedStopwords == null) {
        	BufferedReader br = null;
    		try {
    			br = new BufferedReader(new FileReader(JGraphTTest.class.getResource("/stopwords.txt").getPath().substring(1)));
    		} catch (FileNotFoundException e) {
    			e.printStackTrace();
    		}
        	try {
        	    StringBuilder sb = new StringBuilder();
        	    String line = br.readLine();

        	    while (line != null) {
        	        sb.append(line);
        	        sb.append(stopWordSeparator);
        	        line = br.readLine();
        	    }
        	    linkedStopwords = sb.toString();
        	} catch (IOException e) {
    			e.printStackTrace();
    		} finally {
        	    try {
    				br.close();
    			} catch (IOException e) {
    				e.printStackTrace();
    			}
        	}
    	}
	}

	/**
     * {@inheritDoc}
     */
    @Override
    public void init()
    {
//		String paragraph = "The quick brown fox jumps over the lazy dog. "
//				+ "The quick brown fox jumps over the lazy cat.";
    	
		String paragraph = "The wife of a former U.S. president Bill Clinton Hillary Clinton visited China last Monday. "
				+ "Hillary Clinton wanted to visit China last month but postponed her plans till Monday last week. "
				+ "Hillary Clinton paid a visit to the People Republic of China on Monday. "
				+ "Last week the Secretary of State Ms. Clinton visited Chinese officials.";
		
		// creates a StanfordCoreNLP object, with POS tagging, lemmatization,
		// NER, parsing, and coreference resolution
		Properties props = new Properties();
		/**
		 * tokenize: divides text into a sequence of words
		 * ssplit: Splits a sequence of tokens into sentences(?)
		 * pos: Labels tokens with their POS tag
		 * lemma: distinguish words in different forms, like plural form and tenses
		 * ner: Named Entity Recognizer
		 * parse: syntax analysis
		 * dcoref: coreference resolution
		 */
//		props.put("annotators", "tokenize, ssplit, pos, lemma, ner, parse, dcoref");
		// ssplit depends on tokenize
		props.put("annotators", "tokenize, ssplit, pos");
		StanfordCoreNLP pipeline = new StanfordCoreNLP(props);

		// create an empty Annotation just with the given text
		Annotation document = new Annotation(paragraph);

		// run all Annotators on this text
		pipeline.annotate(document);

		// these are all the sentences in this document
		// a CoreMap is essentially a Map that uses class objects as keys and
		// has values with custom types
		List<CoreMap> sentences = document.get(SentencesAnnotation.class);

		List<ArrayList<WordPosTuple>> sentenceList = new ArrayList<>();
		
		ArrayList<WordPosTuple> startTuple = new ArrayList<>();
		startWord = new WordPosTuple("-start-", "-start-");
		startTuple.add(startWord);
    	wordGraph.addVertex(startWord);
		ArrayList<WordPosTuple> endTuple = new ArrayList<>();
		endWord = new WordPosTuple("-end-", "-start-");
		endTuple.add(endWord);
    	wordGraph.addVertex(endWord);
		
		// add start word
		sentenceList.add(startTuple);
		
		for (CoreMap sentence : sentences) {
			ArrayList<WordPosTuple> wordPosTupleList = new ArrayList<>();
			// traversing the words in the current sentence
			// a CoreLabel is a CoreMap with additional token-specific methods
			for (CoreLabel token : sentence.get(TokensAnnotation.class)) {
				// this is the text of the token
				String word = token.get(TextAnnotation.class);
				// this is the POS tag of the token
				String pos = token.get(PartOfSpeechAnnotation.class);
				// this is the NER label of the token
//				String ne = token.get(NamedEntityTagAnnotation.class);
				if (StringUtils.isPunct(word)) {
					wordPosTupleList.add(new WordPosTuple(word, punct_tag));
				} else {
					wordPosTupleList.add(new WordPosTuple(word, pos));
				}
			}
			System.out.println(wordPosTupleList);
			sentenceList.add(wordPosTupleList);
		}
		
		// add end word
		sentenceList.add(endTuple);
        
		term_freq = compute_statistics(sentenceList);
		
		build_graph(sentenceList);
		
        // note directed edges are printed as: (<v1>,<v2>)
        System.out.println(wordGraph);
        
        // create a JGraphT graph
        ListenableGraph<WordPosTuple, CustomDirectedWeightedEdge> g = wordGraph;

        // create a visualization using JGraph, via an adapter
        jgAdapter = new JGraphModelAdapter<>(g);

        JGraph jgraph = new JGraph(jgAdapter);

        adjustDisplaySettings(jgraph);
        getContentPane().add(jgraph);
        resize(DEFAULT_SIZE);
        
		Set<CustomDirectedWeightedEdge> allEdges = wordGraph.edgeSet();
		for (DefaultWeightedEdge edge : allEdges) {
			System.out.println(edge.toString());
		}
    }

	private void build_graph(List<ArrayList<WordPosTuple>> sentenceList) {
		for (int i = 0; i<sentenceList.size(); i++) {
			ArrayList<WordPosTuple> curr_sentence = sentenceList.get(i);
			
			int sentence_len = sentenceList.get(i).size();
			
			WordPosTuple[] mapping = new WordPosTuple[sentence_len];
			/**
			 * 1. non-stopwords for which no candidate exists in the graph or
			 * for which an unambiguous mapping is possible or which occur
			 * more than once in the sentence.
			 */
			System.out.println("Step 1...");
			for (int j = 0; j<sentence_len; j++) {
				WordPosTuple currTuple = curr_sentence.get(j);
				String word = currTuple.getText();
				if (isStopword(word) || StringUtils.isPunct(word)) {
					continue;
				}
				
				int k = ambiguous_nodes(currTuple);
				
				// If there is no node in the graph, create one with id = 0
				if (k == 0) {
					currTuple.setId(0);
					currTuple.addInfo(i, j);
					wordGraph.addVertex(currTuple);
					mapping[j] = currTuple;
				} else if (k == 1) {
					// Get the sentences id of this node
					Set<Integer> ids = currTuple.getSentenceIdSet();
					// Update the node in the graph if not same sentence
					if (!ids.contains(i)) {
						currTuple.setId(0);
						currTuple.addInfo(i, j);
						wordGraph.addVertex(currTuple);
						mapping[j] = currTuple;
					// Else Create new node for redundant word
					} else {
						currTuple.setId(1);
						currTuple.addInfo(i, j);
						wordGraph.addVertex(currTuple);
						mapping[j] = currTuple;
					}
				}
			}
			
			/**
			 * 2. non-stopwords for which there are either several possible
			 * candidates in the graph.
			 */
			System.out.println("Step 2...");
			for (int j = 0; j<sentence_len; j++) {
				WordPosTuple currTuple = curr_sentence.get(j);
				String word = currTuple.getText();
				if (isStopword(word) || StringUtils.isPunct(word)) {
					continue;
				}
				
				// If word is not already mapped to a node
				if (mapping[j] == null) {
					WordPosTuple prev_tuple = getPrevious(curr_sentence, j);
					WordPosTuple next_tuple = getNext(curr_sentence, j);
				
					int k = ambiguous_nodes(currTuple);
					
					HashMap<WordPosTuple, Integer> ambinode_overlap = new HashMap<>();
					HashMap<WordPosTuple, Integer> ambinode_frequency = new HashMap<>();
						
					for (int c=0; c<k; c++) {
						WordPosTuple tupleToCheck = new WordPosTuple(currTuple, c);
						Set<WordPosTuple> allTuples = wordGraph.vertexSet();
						for (WordPosTuple t : allTuples) {
							if (t.equals(tupleToCheck)) {
								tupleToCheck = t;
								break;
							}
						}
						Set<WordPosTuple> l_context = get_directed_context(sentenceList, tupleToCheck, "left", false);
						Set<WordPosTuple> r_context = get_directed_context(sentenceList, tupleToCheck, "right", false);
						
						int l_context_count = 0;
						int r_context_count = 0;
						for (WordPosTuple t : l_context) {
							if (t.hasSameTxtAndPos(prev_tuple))
								l_context_count++;
						}
						for (WordPosTuple t : r_context) {
							if (t.hasSameTxtAndPos(next_tuple))
								r_context_count++;
						}
						
						int val = l_context_count+r_context_count;
						
						ambinode_overlap.put(tupleToCheck, val);
						
						ambinode_frequency.put(tupleToCheck, tupleToCheck.getInfoList().size());
					}
						
					boolean found = false;
					WordPosTuple selected = null;
					while (!found) {
						// Select the ambiguous node
						selected = max_index(ambinode_overlap);
						if (ambinode_overlap.get(selected) != null) {
							selected = max_index(ambinode_frequency);
						}
						
						// Get the sentences id of this node
						Set<Integer> ids = new HashSet<>();
						if (selected != null) {
							selected.getSentenceIdSet();
						}
						
						// Test if there is no loop
						Integer integer = (Integer) i;
	                    if (!ids.contains(integer)) {
	                    	found = true;
	                    	break;
	                	// Remove the candidate from the lists
	                    } else {
	                    	ambinode_overlap.remove(selected);
	                    	ambinode_frequency.remove(selected);
	                    }
	                    
	                    // Avoid endless loops
	                    if (ambinode_overlap.isEmpty()) {
	                    	break;
	                    }
					}
			                
	                // Update the node in the graph if not same sentence
					if (found && selected != null) {
						selected.addInfo(i, j);
						wordGraph.addVertex(selected);
						mapping[j] = selected;
					} else {
						// Else create new node for redundant word
						currTuple.setId(k);
						currTuple.addInfo(i, j);
						wordGraph.addVertex(currTuple);
						mapping[j] = currTuple;
					}
				}
			}
			
			/**
			 * 3. map the stopwords to the nodes
			 */
			System.out.println("Step 3...");
			for (int j = 0; j<sentence_len; j++) {
				WordPosTuple currTuple = curr_sentence.get(j);
				String word = currTuple.getText();
				
				// If *NOT* stopword, continues
				if (!isStopword(word)) {
					continue;
				}
				
				// Find the number of ambiguous nodes in the graph
				int k = ambiguous_nodes(currTuple);
				
				// If there is no node in the graph, create one with id = 0
				if (k == 0) {
					// Add the node in the graph
					currTuple.setId(0);
					currTuple.addInfo(i, j);
					wordGraph.addVertex(currTuple);
					mapping[j] = currTuple;
				// Else find the node with overlap in context or create one
				} else {
					WordPosTuple prev_tuple = getPrevious(curr_sentence, j);
					WordPosTuple next_tuple = getNext(curr_sentence, j);
					
					HashMap<WordPosTuple, Integer> ambinode_overlap = new HashMap<>();
					
					for (int c=0; c<k; c++) {
						WordPosTuple tupleToCheck = new WordPosTuple(currTuple, c);
						Set<WordPosTuple> allTuples = wordGraph.vertexSet();
						for (WordPosTuple t : allTuples) {
							if (t.equals(tupleToCheck)) {
								tupleToCheck = t;
								break;
							}
						}
						Set<WordPosTuple> l_context = get_directed_context(sentenceList, tupleToCheck, "left", true);
						Set<WordPosTuple> r_context = get_directed_context(sentenceList, tupleToCheck, "right", true);
						
						int l_context_count = 0;
						int r_context_count = 0;
						for (WordPosTuple t : l_context) {
							if (t.hasSameTxtAndPos(prev_tuple))
								l_context_count++;
						}
						for (WordPosTuple t : r_context) {
							if (t.hasSameTxtAndPos(next_tuple))
								r_context_count++;
						}
						
						int val = l_context_count+r_context_count;
						
						ambinode_overlap.put(tupleToCheck, val);
					}
					
					// Get best overlap candidate
					WordPosTuple selected = max_index(ambinode_overlap);
					
					// Get the sentences id of this node
					Set<Integer> ids = new HashSet<>();
					if (selected != null) {
						ids = selected.getSentenceIdSet();
					}
					
					// Update the node in the graph if not same sentence and
					// there is at least one overlap in context
					Integer integer = (Integer) i;
					if (selected != null && !ids.contains(integer) && !ambinode_overlap.isEmpty()) {
						// Update the node in the graph
						selected.addInfo(i, j);
						mapping[j] = selected;
					// Else create a new node
					} else {
						// Add the node in the graph
						currTuple.setId(k);
						currTuple.addInfo(i, j);
						wordGraph.addVertex(currTuple);
						mapping[j] = currTuple;
					}
				}
			}
			
			/**
			 * 4. lastly map the punctuation marks to the nodes
			 */
			System.out.println("Step 4...");
			for (int j = 0; j<sentence_len; j++) {
				WordPosTuple currTuple = curr_sentence.get(j);
				
				if (!StringUtils.isPunct(currTuple.getText())) {
					continue;
				}
				
				int k = ambiguous_nodes(currTuple);
				
				if (k == 0) {
					currTuple.setId(0);
					currTuple.addInfo(i, j);
					wordGraph.addVertex(currTuple);
					mapping[j] = currTuple;
				} else {
					WordPosTuple prevTuple = getPrevious(curr_sentence, j);
					WordPosTuple nextTuple = getNext(curr_sentence, j);
					
					HashMap<WordPosTuple, Integer> ambinode_overlap = new HashMap<>();
					
					for (int c=0; c<k; c++) {
						WordPosTuple tupleToCheck = new WordPosTuple(currTuple, c);
						Set<WordPosTuple> allTuples = wordGraph.vertexSet();
						for (WordPosTuple t : allTuples) {
							if (t.equals(tupleToCheck)) {
								tupleToCheck = t;
								break;
							}
						}
						Set<WordPosTuple> l_context = get_directed_context(sentenceList, tupleToCheck, "left", true);
						Set<WordPosTuple> r_context = get_directed_context(sentenceList, tupleToCheck, "right", true);
						
						int l_context_count = 0;
						int r_context_count = 0;
						for (WordPosTuple t : l_context) {
							if (t.hasSameTxtAndPos(prevTuple))
								l_context_count++;
						}
						for (WordPosTuple t : r_context) {
							if (t.hasSameTxtAndPos(nextTuple))
								r_context_count++;
						}
						
						int val = l_context_count+r_context_count;
						
						ambinode_overlap.put(tupleToCheck, val);
					}
					
					WordPosTuple selected = max_index(ambinode_overlap);
					
					Set<Integer> ids = new HashSet<>();
					if (selected != null) {
						selected.getSentenceIdSet();
					}

					Integer integer = (Integer) i;
					if (selected != null && !ids.contains(integer) && ambinode_overlap.size() > 1) {
						selected.addInfo(i, j);
						wordGraph.addVertex(selected);
						mapping[j] = selected;
					} else {
						currTuple.setId(k);
						currTuple.addInfo(i, j);
						wordGraph.addVertex(currTuple);
						mapping[j] = currTuple;
					}
				}
			}
			
			/**
			 * 4. Connects the mapped words with directed edges
			 */
			for (int j = 1; j<mapping.length; j++) {
				System.out.printf("(%s: %d, %s: %d)\n",mapping[j-1],mapping[j-1].getInfoList().size(), mapping[j], mapping[j].getInfoList().size());
				wordGraph.addEdge(mapping[j-1], mapping[j]);
			}
		}
		
		// Assigns a weight to each node in the graph
		Set<CustomDirectedWeightedEdge> allEdges = wordGraph.edgeSet();
		for (CustomDirectedWeightedEdge edge : allEdges) {
			// Get the list of (sentence_id, pos_in_sentence) for node1
			WordPosTuple source = wordGraph.getEdgeSource(edge);
			// Get the list of (sentence_id, pos_in_sentence) for node2
			WordPosTuple target = wordGraph.getEdgeTarget(edge);
			
			double edge_weight = get_edge_weight(sentenceList, source, target);
			wordGraph.setEdgeWeight(edge, edge_weight);
		}
	}
	
	private double get_edge_weight(List<ArrayList<WordPosTuple>> sentenceList, WordPosTuple source, WordPosTuple target) {
		ArrayList<InfoTuple> info1 = source.getInfoList();
		ArrayList<InfoTuple> info2 = target.getInfoList();
		
		int freq1 = info1.size();
		int freq2 = info2.size();
		
		ArrayList<Double> diff = new ArrayList<>(); 
		
		for (int s=0; s<sentenceList.size(); s++) {
			ArrayList<Integer> pos_i_in_s = new ArrayList<>();
			ArrayList<Integer> pos_j_in_s = new ArrayList<>();
			
			for (InfoTuple itp : info1) {
				if (itp.getSentenceId() == s) {
					pos_i_in_s.add(itp.getPosInSentence());
				}
			}
			
			for (InfoTuple itp : info2) {
				if (itp.getSentenceId() == s) {
					pos_j_in_s.add(itp.getPosInSentence());
				}
			}
			
			ArrayList<Double> all_diff_pos_i_j = new ArrayList<>();
			
			// Loop over all the i, j couples
			for (int x=0; x<pos_i_in_s.size(); x++) {
				for (int y=0; y<pos_j_in_s.size(); y++) {
					double diff_i_j = pos_i_in_s.get(x) - pos_j_in_s.get(y);
					// Test if word i appears *BEFORE* word j in s
					if (diff_i_j<0) {
						all_diff_pos_i_j.add(-1.0*diff_i_j);
					}
				}
			}
			
			// Add the minimum distance to diff (i.e. in case of multiple
			// occurrences of i or/and j in sentence s), 0 otherwise.
			if (all_diff_pos_i_j.size() > 0) {
				diff.add(1.0/Collections.min(all_diff_pos_i_j));
			} else {
				diff.add(0.0);
			}
		}
		
		int weight1 = freq1;
		int weight2 = freq2;
		
		return ( (freq1 + freq2) / diff.stream().mapToDouble(Double::doubleValue).sum() ) / (weight1 * weight2);
	}

	private static Set<WordPosTuple> get_directed_context(List<ArrayList<WordPosTuple>> sentenceList, WordPosTuple wpt, String dir, boolean non_pos) {
		// Define the context containers
		Set<WordPosTuple> l_context = new HashSet<>();
		Set<WordPosTuple> r_context = new HashSet<>();
		
		for (InfoTuple itp : wpt.getInfoList()) {
			// For all the sentence/position tuples
			WordPosTuple prev = getPrevious(sentenceList.get(itp.getSentenceId()), itp.getPosInSentence());
			WordPosTuple next = getNext(sentenceList.get(itp.getSentenceId()), itp.getPosInSentence());
			
			if (non_pos) {
				if (prev != null && !isStopword(prev.getText())) {
					l_context.add(prev);
				}				
				if (next != null && !isStopword(next.getText())) {
					r_context.add(next);
				}
			}
			
			// Returns the left (previous) context
			if (dir.equals("left")) {
				return l_context;
			// Returns the right (next) context
			} else if (dir.equals("right")) {
				return r_context;
			// Returns the whole context
			} else {
				l_context.addAll(r_context);
				return l_context;
			}
		}
		
		return l_context;
	}

	private static WordPosTuple max_index(HashMap<WordPosTuple, Integer> map) {
		if (map.isEmpty()) return null;
		
		int max_val = 0;
		WordPosTuple max_tuple = null;
   		for (Entry<WordPosTuple, Integer> entry : map.entrySet()) {
   			if (entry.getValue() > max_val) {
   				max_val = entry.getValue();
   				max_tuple = entry.getKey();
   			}
   		}
		return max_tuple;
	}

	private int ambiguous_nodes(WordPosTuple currTuple) {
		int nb_node = 0;
		Set<WordPosTuple> allTuples = wordGraph.vertexSet();
		WordPosTuple tupleToCheck = currTuple;
		while (allTuples.contains(tupleToCheck)) {
			nb_node += 1;
			tupleToCheck = new WordPosTuple(tupleToCheck, nb_node);
		}
		return nb_node;
	}

	private static HashMap<WordPosTuple, Integer> compute_statistics(List<ArrayList<WordPosTuple>> sentenceList) {
		// key: tuple of word and pos
		// value: sentences which contain the tuple
		HashMap<WordPosTuple, List<ArrayList<WordPosTuple>>> terms = new HashMap<>();
		
		// loop over sentences
		for (ArrayList<WordPosTuple> sentence : sentenceList) {
			// for each tuple
			for (WordPosTuple tuple : sentence) {
				if (!terms.containsKey(tuple)) {
					// create new entry for new tuple
					List<ArrayList<WordPosTuple>> value = new ArrayList<>();
					value.add(sentence);
					terms.put(tuple, value);
				} else {
					// add sentence entry for existing tuple
					List<ArrayList<WordPosTuple>> value = terms.get(tuple);
					value.add(sentence);
					terms.put(tuple, value);
				}
			}
		}
		
		HashMap<WordPosTuple, Integer> term_freq = new HashMap<>();
		// for each tuple
   		for (Map.Entry<WordPosTuple, List<ArrayList<WordPosTuple>>> entry : terms.entrySet()) {
   			// compute frequency
   			WordPosTuple tuple = entry.getKey();
   			List<ArrayList<WordPosTuple>> sentences = entry.getValue();
   			term_freq.put(tuple, sentences.size());
   		}
   		return term_freq;
	}

	private void adjustDisplaySettings(JGraph jg)
    {
        jg.setPreferredSize(DEFAULT_SIZE);

        Color c = DEFAULT_BG_COLOR;
        String colorStr = null;

        try {
            colorStr = getParameter("bgcolor");
        } catch (Exception e) {
        }

        if (colorStr != null) {
            c = Color.decode(colorStr);
        }

        jg.setBackground(c);
    }

    @SuppressWarnings("unchecked")
    private void positionVertexAt(Object vertex, int x, int y)
    {
        DefaultGraphCell cell = jgAdapter.getVertexCell(vertex);
        AttributeMap attr = cell.getAttributes();
        Rectangle2D bounds = GraphConstants.getBounds(attr);

        Rectangle2D newBounds = new Rectangle2D.Double(x, y, bounds.getWidth(), bounds.getHeight());

        GraphConstants.setBounds(attr, newBounds);

        AttributeMap cellAttr = new AttributeMap();
        cellAttr.put(cell, attr);
        jgAdapter.edit(cellAttr, null, null, null);
    }

    /**
     * a listenable directed multigraph that allows loops and parallel edges.
     */
    private static class ListenableDirectedWeightedGraph<V, E>
        extends DefaultListenableGraph<V, E>
        implements DirectedGraph<V, E>, WeightedGraph<V, E>
    {
        private static final long serialVersionUID = 1L;

        ListenableDirectedWeightedGraph(Class<E> edgeClass)
        {
            super(new DefaultDirectedWeightedGraph<>(edgeClass));
        }
    }
	
	public static WordPosTuple getNext(ArrayList<WordPosTuple> wordList, int idx) {
//	    int idx = wordList.indexOf(currWord);
	    if (idx < 0 || idx+1 == wordList.size()) return null;
	    return wordList.get(idx + 1);
	}

	public static WordPosTuple getPrevious(ArrayList<WordPosTuple> wordList, int idx) {
//	    int idx = wordList.indexOf(currWord);
	    if (idx <= 0) return null;
	    return wordList.get(idx - 1);
	}

	private static boolean isStopword(String word) {
		ArrayList<String> stopwordList = new ArrayList<>(Arrays.asList(linkedStopwords.split(",")));
		if (stopwordList.contains(word)) {
			return true;
		}
		return false;
	}
}
