package ger.maddison;

import java.util.List;

public interface Node {
	boolean isTerminal();
	
	List<Node> getChildren();
	String getName();
}
