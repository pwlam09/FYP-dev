package com.main;

public class InfoTuple {
	private int sentence_id;
	private int posInSentence;
	
	public InfoTuple(int sentence_id, int posInSentence) {
		this.sentence_id = sentence_id;
		this.posInSentence = posInSentence;
	}

	public int getSentenceId() {
		return sentence_id;
	}

	public int getPosInSentence() {
		return posInSentence;
	}
}
