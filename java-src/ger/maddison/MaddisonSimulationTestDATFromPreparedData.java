package ger.maddison;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import ger.maddison.MaddisonSimulationTestDAT.LabelingGandAndLoss;
import ger.tree.DoubleArrayTree;

public class MaddisonSimulationTestDATFromPreparedData {
	public static int testCall0() {
		return 1;
	}

	public static int testCall1() {
		System.out.println("Test message");
		System.err.println("Test error message");
		return 10;
	}

//	public static void main(String[] args) throws IOException {
//		// File dir = new File("precalculatedTreeLabelings");
//		File dir = new File("./generatedSamples_2022_03_15/");
//		load(69, 4, dir);
//		// load(0, 73, dir);
//	}

	private static final Pattern regexp1 = Pattern.compile("NN(\\d+)-[a-z]+\\.ser\\.zip");

	static ArrayList<LabelingGandAndLoss> load(int nGain, int nLoss, File dir) throws IOException {
//	System.out.println("DEBUG start loading 2");
		// Pattern regexp1 =
		// Pattern.compile("NN(\\d+)_GAINS(\\d+)_N(\\d+)-[a-z]+\\.ser");
//	System.out.println("DEBUG start loading 2.1");

//	System.out.println("DEBUG start loading 3");

		ArrayList<LabelingGandAndLoss> result = new ArrayList<LabelingGandAndLoss>();

		int nn = nGain + nLoss;
		String[] fileList = dir.list();

//	System.out.println("DEBUG start loading 4");

		int sampleCnt1 = 0;
		int sampleCnt0 = 0;

//	System.out.println("DEBUG filelist.length=" + fileList.length);
//	System.out.flush();
		for (int i = 0; i < fileList.length; i++) {
			Matcher m = regexp1.matcher(fileList[i]);
			if (m.matches()) {
				int nn1 = Integer.parseInt(m.group(1));
				// int nGain1 = Integer.parseInt(m.group(2));
				// int n1 = Integer.parseInt(m.group(3));

				if (nn1 == nn) {
					System.out.println(fileList[i]);

					File f = new File(dir, fileList[i]);
					try {
						Map<Integer, ArrayList<LabelingGandAndLoss>> map = PrepareMaddisonSimulation
								.loadFromZipFile(f);
						System.out.println(map.keySet());
						// StringBuilder sb = new StringBuilder();
						// sb.append("[ ");
						// for (Map.Entry<Integer, ArrayList<LabelingGandAndLoss>> entry :
						// map.entrySet())
						// {
						// sb.append(entry.getKey()).append("(").append(entry.getValue().size()).append("),
						// ");
						// }
						// sb.setLength(sb.length() - 2);
						// sb.append("]");
						// System.out.println(sb);

						ArrayList<LabelingGandAndLoss> list1 = map.get(nGain);
						if (list1 != null) {
							result.addAll(list1);
							sampleCnt0 += list1.size();
						}
						ArrayList<LabelingGandAndLoss> list2 = map.get(nLoss);
						if (list2 != null) {
							reverseAll(list2);
							result.addAll(list2);
							sampleCnt1 += list2.size();
						}

					} catch (java.io.IOException e) {
						System.out.println("HIBAS:" + fileList[i]);
						e.printStackTrace();
					}

				}

			} else {
				System.out.println("unknown: '" + fileList[i] + "'");
			}

		}
		System.out.println("Size of random sample set:" + result.size() + " = " + sampleCnt0 + "+" + sampleCnt1);

		return result;
	}

	public static Precalculated precalculate(int nGain, int nLoss, DoubleArrayTree tree, File dir) throws IOException {
//	System.out.println("DEBUG factory function called");
		try {
			return new Precalculated(nGain, nLoss, tree, dir);
		} catch (IOException e) {
			e.printStackTrace();
			throw (e);
		} catch (RuntimeException e) {
			e.printStackTrace();
			throw (e);
		}
	}

	public static class Precalculated {

		final ArrayList<LabelingGandAndLoss> labelingList1;
		final DoubleArrayTree tree;
		int nGain;
		int nLoss;

		public Precalculated(int nGain, int nLoss, DoubleArrayTree tree, File dir) throws IOException {

			this.tree = tree;
			this.nGain = nGain;
			this.nLoss = nLoss;

//	    System.out.println("DEBUG start loading");
			ArrayList<LabelingGandAndLoss> labelingList1 = load(nGain, nLoss, dir);
			ArrayList<LabelingGandAndLoss> labelingList2 = load(nLoss, nGain, dir);
			reverseAll(labelingList2);
			labelingList1.addAll(labelingList2);
			labelingList2 = null;
			this.labelingList1 = labelingList1;

//	    System.out.println("DEBUG constructor finished");
		}

		public int size() {
			return labelingList1.size();
		}

		@Override
		public String toString() {
			return "PrecalculatesSamples[gains=" + nGain + ", losses=" + nLoss + "sample_cnt=" + size() + "tree_size="
					+ tree.size() + "]";
		}

		/**
		 * @return  (nGain + 1) * (nLoss + 1) sized array, megszamolva, hogy az egyes (gain, loss) parokhoz hany talalat van a szimulaciobol
		 */
		public int[] calulateTable(boolean[] idxOfBlackNodes) {

			int[] resultTbl = new int[(nGain + 1) * (nLoss + 1)];

			for (LabelingGandAndLoss gl : labelingList1) {
				int cntOfGainsInBlackArea = 0;
				int cntOfLossesInBlackArea = 0;

				for (int x : gl.gainArray) {
					cntOfGainsInBlackArea += idxOfBlackNodes[x] ? 1 : 0;
				}
				for (int x : gl.lossArray) {
					cntOfLossesInBlackArea += idxOfBlackNodes[x] ? 1 : 0;
				}

				resultTbl[cntOfGainsInBlackArea * (1 + nLoss) + cntOfLossesInBlackArea]++;
			}
			return resultTbl;

		}

	}

	public static int[] calulateTable(int nGain, int nLoss, DoubleArrayTree tree, boolean[] idxOfBlackNodes, File dir)
			throws IOException {
		Precalculated p = new Precalculated(nGain, nLoss, tree, dir);
		return p.calulateTable(idxOfBlackNodes);
	}

	static LabelingGandAndLoss reverse(LabelingGandAndLoss lg) {
		// Forditott sorrendbe adom be a parametereket !!!
		return new LabelingGandAndLoss(lg.lossArray, lg.gainArray);
	}

	static void reverseAll(ArrayList<LabelingGandAndLoss> list) {
		final int n = list.size();
		for (int i = 0; i < n; i++) {
			list.set(i, reverse(list.get(i)));
		}
	}

}