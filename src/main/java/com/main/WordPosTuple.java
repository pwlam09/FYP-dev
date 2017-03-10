package com.main;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Set;

public class WordPosTuple {
	private String word;
	private String pos;
	private int id;
	private ArrayList<InfoTuple> infoList;
	
	public WordPosTuple(String word, String pos) {
		this.word = word;
		this.pos = pos;
		this.infoList = new ArrayList<>();
	}

	public WordPosTuple(WordPosTuple tupleToCheck, int id) {
		this.word = tupleToCheck.word;
		this.pos = tupleToCheck.pos;
		this.infoList = new ArrayList<>();
		this.id = id;
	}

	public String getText() {
		return word;
	}

	public int getId() {
		return id;
	}
	
	public boolean hasSameTxtAndPos(WordPosTuple otherWord) {
	    if (word.toLowerCase().equals(otherWord.word.toLowerCase()) && pos.equals(otherWord.pos)) {
	    	return true;
	    } else {
	    	return false;
	    }
	}
	
	@Override
	public boolean equals(Object obj) {
	    if (obj == null) return false;
	    if (obj == this) return true;
	    if (!(obj instanceof WordPosTuple)) return false;
	    
	    WordPosTuple otherWord = (WordPosTuple) obj;
	    if (word.toLowerCase().equals(otherWord.word.toLowerCase()) && pos.equals(otherWord.pos) && id == otherWord.id) {
	    	return true;
	    } else {
	    	return false;
	    }
	}

	@Override
	public int hashCode() {
//        int hash = 7;
//        hash = 47 * hash + Objects.hashCode(this.text);
//        hash = 47 * hash + Objects.hashCode(this.pos);
        return id;
	}

	@Override
	public String toString() {
		return String.format("%s_%s hash: %d", word, pos, id);
	}
	
	public void setId(int id) {
		this.id = id;
	}
	
	public void addInfo(int sentence_id, int position_in_sentence) {
		this.infoList.add(new InfoTuple(sentence_id, position_in_sentence));
	}

	public ArrayList<InfoTuple> getInfoList() {
		return infoList;
	}

	public Set<Integer> getSentenceIdSet() {
		Set<Integer> sidList = new HashSet<>();
		for (InfoTuple t : infoList) {
			sidList.add(t.getSentenceId());
		}
		return sidList;
	}
}
