package com.main;
import java.awt.FlowLayout;
import java.awt.Image;
import java.awt.image.BufferedImage;
import java.awt.image.DataBufferByte;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;

import javax.swing.ImageIcon;
import javax.swing.JFrame;
import javax.swing.JLabel;

import org.opencv.core.Mat;

/**
 * @author PuiWa
 *
 */
public class Main {
	/**
	 * deeplearning4j crash with slib because of different slf4j dependency
	 */
	
	private static String videoPath = "";
	
	private static void main(String[] args) throws Exception {
		/**
		 * python trial
		 */

		int number1 = 10;
		int number2 = 32;
		try {
			String prg = "import sys\nprint int(sys.argv[1])+int(sys.argv[2])\n";
			BufferedWriter out = new BufferedWriter(new FileWriter("test1.py"));
			out.write(prg);
			out.close();

			ProcessBuilder pb = new ProcessBuilder("python","test1.py",""+number1,""+number2);
			Process p = pb.start();

			BufferedReader in = new BufferedReader(new InputStreamReader(p.getInputStream()));
			int ret = new Integer(in.readLine()).intValue();
			System.out.println("value is : " + ret);
		} catch (Exception e) {
			System.out.println(e);
		}
		
		ProcessBuilder pb = new ProcessBuilder("python","test1.py",""+number1,""+number2);
		Process p = pb.start();
		BufferedReader in = new BufferedReader(new InputStreamReader(p.getInputStream()));
		int ret = new Integer(in.readLine()).intValue();
		System.out.println("value is : " + ret);
		/**
		 * slib attempt
		 */
//		// Location of WordNet Data
//        String dataloc = "C:/PP_program/wordnet/dict/";
//
//        // We create the graph
//        URIFactory factory = URIFactoryMemory.getSingleton();
//        URI guri = factory.getURI("http://graph/wordnet/");
//        G wordnet = new GraphMemory(guri);
//
//        // We load the data into the graph
//        GraphLoader_Wordnet loader = new GraphLoader_Wordnet();
//
//        GDataConf dataNoun = new GDataConf(GFormat.WORDNET_DATA, dataloc + "data.noun");
//        GDataConf dataVerb = new GDataConf(GFormat.WORDNET_DATA, dataloc + "data.verb");
//        GDataConf dataAdj = new GDataConf(GFormat.WORDNET_DATA, dataloc + "data.adj");
//        GDataConf dataAdv = new GDataConf(GFormat.WORDNET_DATA, dataloc + "data.adv");
//
//        loader.populate(dataNoun, wordnet);
//        loader.populate(dataVerb, wordnet);
//        loader.populate(dataAdj, wordnet);
//        loader.populate(dataAdv, wordnet);
//
//        // We root the graph which has been loaded (this is optional but may be required to compare synset which do not share common ancestors).
//        GAction addRoot = new GAction(GActionType.REROOTING);
//        GraphActionExecutor.applyAction(addRoot, wordnet);
//
//        // This is optional. It just shows which are the synsets which are not subsumed
//        ValidatorDAG validatorDAG = new ValidatorDAG();
//        Set<URI> roots = validatorDAG.getTaxonomicRoots(wordnet);
//        System.out.println("Roots: " + roots);
//
//        // We create an index to map the nouns to the vertices of the graph
//        // We only build an index for the nouns in this example
//        String data_noun = dataloc + "index.noun";
//
//        IndexerWordNetBasic indexWordnetNoun = new IndexerWordNetBasic(factory, wordnet, data_noun);
//
//        // uncomment if you want to show the index, i.e. nouns and associated URIs (identifiers)
//        for (Map.Entry<String, Set<URI>> entry : indexWordnetNoun.getIndex().entrySet()) {
//            System.out.println(entry.getKey() + "\t" + entry.getValue());
//        }
//
//        // We focus on three specific nouns in this example
//        // - iced_coffee	[http://graph/wordnet/07936780]
//        // - instant_coffee	[http://graph/wordnet/07936903]
//        // - green_tea          [http://graph/wordnet/07951392]
//        // We retrive their identifiers
//        Set<URI> uris_iced_coffee = indexWordnetNoun.get("iced_coffee");
//        Set<URI> uris_instant_coffee = indexWordnetNoun.get("instant_coffee");
//        Set<URI> uris_green_tea = indexWordnetNoun.get("green_tea");
//
//        // Note that multiple URIs (identifiers) can be associated to the same noun
//        // In this example we only consider nouns associated to a single URI so we retrieve their URI
//        URI uri_iced_coffee = uris_iced_coffee.iterator().next();
//        URI uri_instant_coffee = uris_instant_coffee.iterator().next();
//        URI uri_green_tea = uris_green_tea.iterator().next();
//
//        // We configure a pairwise semantic similarity measure, 
//        // i.e., a measure which will be used to assess the similarity 
//        // of two nouns regarding their associated vertices in WordNet
//        ICconf iconf = new IC_Conf_Topo(SMConstants.FLAG_ICI_SECO_2004);
//        SMconf measureConf = new SMconf(SMConstants.FLAG_SIM_PAIRWISE_DAG_NODE_LIN_1998);
//        measureConf.setICconf(iconf);
//
//        // We define the engine used to compute the score of the configured measure
//        // several preprocessing will be made prior to the first computation, e.g. to compute the Information Content (IC)
//        // of the vertices. This may take some few secondes
//        SM_Engine engine = new SM_Engine(wordnet);
//
//        // we compute the semantic similarities
//        double sim_iced_coffee_vs_instant_coffee = engine.compare(measureConf, uri_iced_coffee, uri_instant_coffee);
//        double sim_iced_coffee_vs_green_tea = engine.compare(measureConf, uri_iced_coffee, uri_green_tea);
//
//        // That's it
//        System.out.println("sim(iced_coffee,instant_coffee) = " + sim_iced_coffee_vs_instant_coffee);
//        System.out.println("sim(iced_coffee,green_tea)      = " + sim_iced_coffee_vs_green_tea);
//
//        // Which prints
//        // sim(iced_coffee,instant_coffee) = 0.7573022852697784
//        // sim(iced_coffee,green_tea)      = 0.3833914674618656
//		
        
        /**
         * tf-idf attempt
         */
//	    List<String> doc1 = Arrays.asList("Lorem", "ipsum", "dolor", "ipsum", "sit", "ipsum");
//	    List<String> doc2 = Arrays.asList("Vituperata", "incorrupte", "at", "ipsum", "pro", "quo");
//	    List<String> doc3 = Arrays.asList("Has", "persius", "disputationi", "id", "simul");
//	    List<List<String>> documents = Arrays.asList(doc1, doc2, doc3);
//	 
//	    Main calculator = new Main();
//	    double tfidf = calculator.tfIdf(doc1, documents, "ipsum");
//	    System.out.println("TF-IDF (ipsum) = " + tfidf);
        
        /**
         * test corpus xml parse with apache tika
         */
		// Directory where the files are located
//		File dir = new File("C:/PP_file/cityuFYP/eclipseWorkspace/FypNewModules/src/main/resources/subtitles/OpenSubtitles/en/result");
//
//		// Create a FileFilter that matches ".xml" files
//		FileFilter filter = (File file) -> file.isFile() && file.getName().endsWith(".xml");
//
//		// Get pathnames of matching files.
//		File[] paths = dir.listFiles(filter);
//		
//		for (File file : paths) {
//			System.out.println(file.getName());
//		}
//		
//		String filepath = "C:/PP_file/cityuFYP/eclipseWorkspace/FypNewModules/src/main/resources/subtitles/OpenSubtitles/en/result/2_195899_259553_identity.xml";
//		// detecting the file type
//		BodyContentHandler handler = new BodyContentHandler();
//		Metadata metadata = new Metadata();
//		File tikaFile = new File(filepath);
//		System.out.println(tikaFile.getAbsolutePath());
//		FileInputStream inputstream = new FileInputStream(tikaFile);
//		ParseContext pcontext = new ParseContext();
//
//		// Xml parser
//		XMLParser xmlparser = new XMLParser();
//		xmlparser.parse(inputstream, handler, metadata, pcontext);
////		System.out.println("Contents of the document:" + handler.toString());
//		System.out.println("Metadata of the document:");
//		String[] metadataNames = metadata.names();
//
//		for (String name : metadataNames) {
//			System.out.println(name + ": " + metadata.get(name));
//		}
//		System.out.println(metadataNames.length);
		
		/**
		 * stanford nlp attempt
		 */
	    // Tokenize
//		String paragraph = "The quick brown fox jumps over the lazy dog. The quick brown fox jumps over the lazy dog.";
	    
		// creates a StanfordCoreNLP object, with POS tagging, lemmatization,
		// NER, parsing, and coreference resolution
//		Properties props = new Properties();
//		/**
//		 * tokenize: divides text into a sequence of words
//		 * ssplit: Splits a sequence of tokens into sentences(?)
//		 * pos: Labels tokens with their POS tag
//		 * lemma: distinguish words in different forms, like plural form and tenses
//		 * ner: Named Entity Recognizer
//		 * parse: syntax analysis
//		 * dcoref: coreference resolution
//		 */
////		props.put("annotators", "tokenize, ssplit, pos, lemma, ner, parse, dcoref");
//		// ssplit depends on tokenize
//		props.put("annotators", "tokenize, ssplit");
//		StanfordCoreNLP pipeline = new StanfordCoreNLP(props);
//
//		// create an empty Annotation just with the given text
//		Annotation document = new Annotation(paragraph);
//
//		// run all Annotators on this text
//		pipeline.annotate(document);
//
//		// these are all the sentences in this document
//		// a CoreMap is essentially a Map that uses class objects as keys and
//		// has values with custom types
//		List<CoreMap> sentences = document.get(SentencesAnnotation.class);
//
//		for (CoreMap sentence : sentences) {
//			System.out.println(sentence.toString());
//			// traversing the words in the current sentence
//			// a CoreLabel is a CoreMap with additional token-specific methods
//			for (CoreLabel token : sentence.get(TokensAnnotation.class)) {
//				// this is the text of the token
//				String word = token.get(TextAnnotation.class);
//				// this is the POS tag of the token
//				String pos = token.get(PartOfSpeechAnnotation.class);
//				// this is the NER label of the token
//				String ne = token.get(NamedEntityTagAnnotation.class);
//			}
//		}
//	    
//		File file = new File(Main.class.getResource("/stanford/english-left3words-distsim.tagger").getPath().substring(1));
//		if (file.exists()) {
//			System.out.println("file exists!");
//		}
//		
//		// Initialize the tagger
//		MaxentTagger tagger = new MaxentTagger(Main.class.getResource("/stanford/english-left3words-distsim.tagger").getPath().substring(1));
//		
//		// The sample string
////		String sample = "This is a sample text";
//		
//		String sample = "Oh my god. Oh my god.";
//		
//		// create text file for DocumentPreprocessor
//		File newTxtFile = new File("test.txt");
//		newTxtFile.createNewFile();
//		try{
//		    PrintWriter writer = new PrintWriter("test.txt", "UTF-8");
//		    writer.println(sample);
//		    writer.close();
//		} catch (IOException e) {
//			e.printStackTrace();
//		}
//		
//		System.out.println(newTxtFile.getAbsolutePath());
//		
//		// remove repeated sentence
//		DocumentPreprocessor dp = new DocumentPreprocessor(newTxtFile.getAbsolutePath());
//		ArrayList<String> sentences = new ArrayList<>();
//		for (List<HasWord> sentence : dp) {
//			System.out.println(sentence);
//			sentences.add(SentenceUtils.listToString(sentence));
//			System.out.println(SentenceUtils.listToString(sentence));
//		}
//		
//		boolean[] isRepeatedSentence = new boolean[sentences.size()];
//		
//		for (int i=0; i<sentences.size(); i++) {
//			for (int j=0; j<sentences.size(); j++) {
//				// i-j==-1 to check whether two sentences are consecutively repeated
//				// if so, set the index of 2nd sentence to false for later processing
//				if (i!=j && (i-j==-1) && sentences.get(i).equals(sentences.get(j))) {
//					isRepeatedSentence[j]=true;
//					System.out.printf("sentence %d equals sentence %d\n", i, j);
//				}
//			}
//		}
//		
//		for (int i = 0; i < isRepeatedSentence.length; i++) {
//			if (isRepeatedSentence[i]) {
//				System.out.printf("index: %d, repeatedSentences: %s\n", i, sentences.get(i));
//			}
//		}
//		
//		ArrayList<String> nonRepeatedSentences = new ArrayList<>();
//		// copy non-repeated sentences to new Array
//		for (int i=0; i<sentences.size(); i++) {
//			if (!isRepeatedSentence[i]) {
//				nonRepeatedSentences.add(sentences.get(i));
//			}
//		}
//		
//		System.out.println("nonRepeatedSentences: "+nonRepeatedSentences);
//		
//		// remove repeated words
//		PTBTokenizer<CoreLabel> ptbt = new PTBTokenizer<>(new FileReader(newTxtFile.getAbsolutePath()), new CoreLabelTokenFactory(), "");
//		while (ptbt.hasNext()) {
//			CoreLabel label = ptbt.next();
//			System.out.println(label);
//		}
//		
//		// The tagged string
//		String tagged = tagger.tagString(sample);
//		 
//		// Output the result
//		System.out.println(tagged);
//		
//		System.out.println(parseTaggedWords(tagged));
		
		//*********Face and Mouth detection test********************
//		System.load(Main.class.getResource("/opencv_java310.dll").getPath().substring(1));
////		System.out.println(Main.class.getResource("/opencv_java310.dll").getPath().substring(1));
////		System.out.println(Main.class.getResource("/haarcascade_mcs_mouth.xml").getPath().substring(1));
//		CascadeClassifier faceDetector = new CascadeClassifier(Main.class.getResource("/haarcascade_frontalface_default.xml").getPath().substring(1));
//		CascadeClassifier mouthDetector = new CascadeClassifier(Main.class.getResource("/haarcascade_mcs_mouth.xml").getPath().substring(1));
//		if (faceDetector.empty() || mouthDetector.empty()) {
//			System.out.println("xml empty");
//		}
//        Mat image = Imgcodecs.imread(Main.class.getResource("/test_img3.JPG").getPath().substring(1));
// 
//        MatOfRect faceDetections = new MatOfRect();
//        faceDetector.detectMultiScale(image, faceDetections);
//      
//        MatOfRect mouthDetections = new MatOfRect();
//        mouthDetector.detectMultiScale(image, mouthDetections);
//   
//        System.out.println(String.format("Detected %s faces", faceDetections.toArray().length));
//        
//        System.out.println(String.format("Detected %s mouths", mouthDetections.toArray().length));
// 
//    	for (Rect faceRect : faceDetections.toArray()) {
//        	Imgproc.rectangle(image, new Point(faceRect.x, faceRect.y), new Point(faceRect.x + faceRect.width, faceRect.y + faceRect.height),
//                    new Scalar(255, 0, 0));
//    	}
//        
//
//    	for (Rect faceRect : faceDetections.toArray()) {
//    		Rectangle upperface = new Rectangle(faceRect.x, faceRect.y, faceRect.width, faceRect.height/2);
//    		Rectangle lowerface = new Rectangle(faceRect.x, faceRect.y+faceRect.height/2, faceRect.width, faceRect.height/2);
//    		ArrayList<Rect> possibleMouths = new ArrayList<>();
//    		
//    		// find all possible mouths for each face
//	        for (Rect mouthRect : mouthDetections.toArray()) {
//	        	boolean isPossibleMouth = false;
//        		if (faceRect.contains(new Point(mouthRect.x, mouthRect.y)) && faceRect.contains(new Point(mouthRect.x+mouthRect.width, mouthRect.y))) {
//            		Rectangle mouthRectJava = new Rectangle(mouthRect.x, mouthRect.y, mouthRect.width, mouthRect.height);
//            		Rectangle upperfaceIntersection = upperface.intersection(mouthRectJava);
//            		Rectangle lowerfaceIntersection = lowerface.intersection(mouthRectJava);
//            		if ((upperfaceIntersection.width * upperfaceIntersection.height < mouthRectJava.width * mouthRectJava.height * 0.3) &&
//        				(lowerfaceIntersection.width * lowerfaceIntersection.height > mouthRectJava.width * mouthRectJava.height * 0.7)) {
//            			isPossibleMouth = true;
//            		}
//        		}
//	        	if (isPossibleMouth) {
//	        		possibleMouths.add(mouthRect);
//	            	Imgproc.rectangle(image, new Point(mouthRect.x, mouthRect.y), new Point(mouthRect.x + mouthRect.width, mouthRect.y + mouthRect.height),
//	                        new Scalar(0, 255, 0));
//	        	}
//	        }
//	        
//        	if (!possibleMouths.isEmpty()) {
//    	        Rect selectedRect = possibleMouths.get(0);
//    	        // if more than 1 mouths, select the one closest to center of lower part of face
//    	        if (possibleMouths.size() > 1) {
//    	        	double maxDiff = 10000.00;
//    	        	Point2D lowerfaceCentre = new Point2D.Double(lowerface.x+lowerface.width/2, lowerface.y+lowerface.height/2);
//        	        for (Rect possibleMouthRect : possibleMouths) {
//        	        	Point2D mouthCentre = new Point2D.Double(possibleMouthRect.x+possibleMouthRect.width/2, possibleMouthRect.y+possibleMouthRect.height/2);
//        	        	if (Math.abs(mouthCentre.distance(lowerfaceCentre)) < maxDiff) {
//        	        		maxDiff = Math.abs(mouthCentre.distance(lowerfaceCentre));
//        	        		selectedRect = possibleMouthRect;
//        	        	}
//        	        }
//    	        }
//            	Imgproc.rectangle(image, new Point(selectedRect.x, selectedRect.y), new Point(selectedRect.x + selectedRect.width, selectedRect.y + selectedRect.height),
//                        new Scalar(0, 0, 255));
//        	}
//    	}
// 
//        String filename = "output.png";
//        System.out.println(String.format("Writing %s", filename));
//        Imgcodecs.imwrite(filename, image);
	}
	
	public static BufferedImage Mat2BufferedImage(Mat m) {
	    // Fastest code
	    // output can be assigned either to a BufferedImage or to an Image

	    int type = BufferedImage.TYPE_BYTE_GRAY;
	    if ( m.channels() > 1 ) {
	        type = BufferedImage.TYPE_3BYTE_BGR;
	    }
	    int bufferSize = m.channels()*m.cols()*m.rows();
	    byte [] b = new byte[bufferSize];
	    m.get(0,0,b); // get all the pixels
	    BufferedImage image = new BufferedImage(m.cols(),m.rows(), type);
	    final byte[] targetPixels = ((DataBufferByte) image.getRaster().getDataBuffer()).getData();
	    System.arraycopy(b, 0, targetPixels, 0, b.length);  
	    return image;
	}
	 
	 public static void displayImage(Image img2) {
	    //BufferedImage img=ImageIO.read(new File("/HelloOpenCV/lena.png"));
	    ImageIcon icon=new ImageIcon(img2);
	    JFrame frame=new JFrame();
	    frame.setLayout(new FlowLayout());        
	    frame.setSize(img2.getWidth(null)+50, img2.getHeight(null)+50);     
	    JLabel lbl=new JLabel();
	    lbl.setIcon(icon);
	    frame.add(lbl);
	    frame.setVisible(true);
	    frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	 }
	 
	public double tf(List<String> doc, String term) {
		double result = 0;
		for (String word : doc) {
			if (term.equalsIgnoreCase(word))
				result++;
		}
		return result / doc.size();
	}

	public double idf(List<List<String>> docs, String term) {
		double n = 0;
		for (List<String> doc : docs) {
			for (String word : doc) {
				if (term.equalsIgnoreCase(word)) {
					n++;
					break;
				}
			}
		}
		return Math.log(docs.size() / n);
	}
	
	public double tfIdf(List<String> doc, List<List<String>> docs, String term) {
	    return tf(doc, term) * idf(docs, term);
	}
}
