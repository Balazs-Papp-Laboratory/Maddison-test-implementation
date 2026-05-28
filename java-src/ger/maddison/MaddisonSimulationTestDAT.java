package ger.maddison;

import java.io.Serializable;
import java.util.Arrays;
import java.util.Collections;
import java.util.Random;

import ger.tree.DoubleArrayTree;
import ger.tree.IntegerStack;

public class MaddisonSimulationTestDAT {

	private static void calulateTable(int nGain, int nLoss, DoubleArrayTree tree, int iteration, Random random) {
		int nn = nGain + nLoss;

		int[] resultTbl = new int[(nGain + 1) * (nLoss + 1)];
		int[] resultGainCntHistogram = new int[nn + 1];

		// Set<Node> setOfAllNoeds = NodeUtils.collectAllDescendants(tree);
		// Node[] arrayOfAllNodes = setOfAllNoeds.toArray(new
		// Node[setOfAllNoeds.size()]);

		// ------
		// random black node set osszeallitasa TODO ezt parameterkent kell kapni
		boolean[] idxOfBlackNodes = new boolean[tree.size()]; // default values are false
		{
			int cnt = 80;
			while (cnt > 0) {

				int i = random.nextInt(idxOfBlackNodes.length - 1) + 1; // a root -ot ki kell zárni, ezért van a +1

				if (!idxOfBlackNodes[i]) {
					cnt--;
					idxOfBlackNodes[i] = true;
				}
			}
		}
		// ------

		IntegerStack gray = new IntegerStack(tree.size());
		// ezt kesobb fogom felhasznalni, i, és mindig resetelem es ujrafelhasznalom,
		// hogy sporoljak a memori a foglalasokkal

		for (int ii = 0; ii < iteration; ii++) {

			// Csinalok egy nn elemu random node halmazt, ahol a gain-ek és loss-ok lesznek
			//
			// eloszor megjelolom azokat a poziciokat ahol valamilyen valtozas lesz (loss
			// vagy
			// gain), aztán ezen poziciok nodejait bedobom egy halmazba

			boolean[] idxOfChanges = new boolean[tree.size()]; // default values are false
			int cnt = nn;
			while (cnt > 0) {

				int i = random.nextInt(idxOfChanges.length - 1) + 1; // a root -ot ki kell zárni, ezért van a +1

				if (!idxOfChanges[i]) {
					cnt--;
					idxOfChanges[i] = true;
				}
			}

			///////////////////////////////////////////
			// Eza blok bejárja a fát, és megszámolja hány gain és hány loss esemény van
			// akkor, ha a setOfChaneges-ben lévő változások vannak, és a
			// gyökérben 1 a gén értéke.
			int cntOfGains = 0;
			int cntOfLosses = 0; // normális, h ezt nem használom, mert a cntOfGains+cntOfLosses=nn
			// osszefuggesbol kijon
			int cntOfGainsInBlackArea = 0;
			int cntOfLossesInBlackArea = 0;

			{
				gray.clean();
				gray.push(tree.getRoot());
				int[] geneValueMap = new int[tree.size()];
				Arrays.fill(geneValueMap, -1);// TODO ez lehet, h felesleges
				// LinkedHashMap<Node, Node> parentMap = new LinkedHashMap<Node, Node>();
				NwkNode sentinelNode = new NwkNode(Collections.<Node>emptyList(), "DUMMY_NODE", 0.0); // ez trukkozes: berakkok
				// egy ál-node-ot a Root
				// szulojenek
				// parentMap.put(tree, sentinelNode);
				geneValueMap[tree.getRoot()] = 1; // feltesszuk, hogy a fa gyokereben 1 a gen értéke TODO azt hiszem ez
				// magatol beall

				while (!gray.isEmpty()) {
					int x = gray.pop();
					int parent = tree.getParent(x);
					int parentValue = parent == -1 ? 1 : geneValueMap[parent];
					boolean isChange = idxOfChanges[x];
					Integer value;
					if (isChange) {
						int isBlack = idxOfBlackNodes[x] ? 1 : 0;
						if (parentValue == 0) {
							value = 1;
							cntOfGains++;
							cntOfGainsInBlackArea += isBlack;
						} else {
							value = 0;
							cntOfLosses++;
							cntOfLossesInBlackArea += isBlack;
						}

					} else {
						value = parentValue;
					}

					// boolean notAlreadyContaind = (geneValueMap[x]== -1);
					geneValueMap[x] = value;

					gray.pushChildren(tree, x); // itt mondom, ha fa obbi reszet is nezze meg

				} // while
			} // fabejárós blokk

			resultGainCntHistogram[cntOfGains]++;

			if (cntOfGains == nGain) {
				// ez a normális ág
				resultTbl[cntOfGainsInBlackArea * (1 + nLoss) + cntOfLossesInBlackArea]++;
			}

			if (cntOfGains == nLoss) {
				// ebben az ágban a gain és loss szerepe felcserélődik
				resultTbl[cntOfLossesInBlackArea * (1 + nLoss) + cntOfGainsInBlackArea]++;
			}

			// System.out.println(cntOfLosses +" "+ cntOfGains);
		}

		System.out.println(Arrays.toString(resultGainCntHistogram));
		System.out.println(resultGainCntHistogram[nGain]);
		System.out.println(((double) iteration) / resultGainCntHistogram[nGain]);
		// TODO kéne valamit csinálni, h mindkét ág ne fusson le, ha nGain==nLoss
	}

	static class LabelingGandAndLoss implements Serializable {
		private static final long serialVersionUID = 3377161367187991572L;

		public final int[] gainArray;
		public final int[] lossArray;

		public LabelingGandAndLoss(int[] gainArray, int[] lossArray) {
			super();
			this.gainArray = gainArray;
			this.lossArray = lossArray;
		}

		@Override
		public String toString() {
			return "GAINS" + Arrays.toString(gainArray) + " LOSSES" + Arrays.toString(lossArray) + ")";
		}
	}

}