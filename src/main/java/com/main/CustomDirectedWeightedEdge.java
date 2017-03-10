package com.main;

import org.jgrapht.graph.DefaultWeightedEdge;

public class CustomDirectedWeightedEdge extends DefaultWeightedEdge {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	@Override
	public String toString() {
		return String.format("(%s->%s, w: %d)", getSource().toString(), getTarget().toString(), (int)Math.floor(getWeight()));
	}
}
