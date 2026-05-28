package ger.maddison;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.util.Arrays;
import java.util.Random;

import ger.maddison.MaddisonSimulationTestDATFromPreparedData.Precalculated;
import ger.tree.DATUtils;
import ger.tree.DoubleArrayTree;

public class DoMaddisonTestForAllPairs {
	static final int ABSENT = 0;
	static final int GAIN = 1;
	static final int LOSS = 2;
	static final int PRESENT = 3;

	private static final String HEADER_PART_OF_P_VALUES = "expedted_number_of_black_gains,expedted_number_of_black_loss,estimated_median_of_black_gains,estimated_median_of_black_loss,"
			+ "p_value_less_gains_then_expected,p_value_more_gains_then_expected,p_value_less_loss_then_expected,p_value_more_loss_then_expected,";
	private static final String HEADER = "cluster1,cluster2,black_area_ratio,sample_size_of_simulation," + //
			"G1_absent,G1_gain,G1_loss,G1_present," + //
			"actual_black_absent,actual_black_gain,actual_black_loss,actual_black_present," + //
			"CM_G1_absent_G2_absent,CM_G1_gain_G2_absent,CM_G1_loss_G2_absent,CM_G1_present_G2_absent,CM_G1_absent_G2_gain,CM_G1_gain_G2_gain,CM_G1_loss_G2_gain,CM_G1_present_G2_gain,CM_G1_absent_G2_loss,CM_G1_gain_G2_loss,CM_G1_loss_G2_loss,CM_G1_present_G2_loss,CM_G1_absent_G2_present,CM_G1_gain_G2_present,CM_G1_loss_G2_present,CM_G1_present_G2_present,"
			+ //
			HEADER_PART_OF_P_VALUES + //
			"table_undersore_separated" + //
			"\n";

	public static void main(String[] args) throws IOException {

		System.out.println(Arrays.toString(args));
		if (args.length != 6 & args.length != 7) {
			System.out.println("7 parameter van,\n" //
					+ "   [1] a hasznát magok száma\n"//
					+ "   [2] fát tartalmazo nwk file\n"//
					+ "   [3] A gyoker elem neve, ha egy reszfat akarunk. Ha nem akkor adj 'root' erteket\n"
					+ "   [4] A 0/1 matrix fileneve, pl data/gene_presence_mx-big-DATorder.tsv\n"
					+ "   [5] Az eloregeneralt random mintak helye, pl \"./linked_dirs/generatedSamples\"\n"
					+ "   [6] A kimeneti könyvtar, ahova az eredmeny tsv fileokat rakja majd., pl \"data-processed/Maddison-resut-java\"\n"
					+ "   [7] a fekete nodeok kiválaszása, pl 11, vagy a default erteke: 01+11+10\n");
		} else {
			final int N = Integer.parseInt(args[0]);// paralell magok száma
			String treeNwkFilename = args[1].trim();
			String rootNodeName = args[2].trim();
			final String filepathOfMatrix = args[3].trim();
			// final String filepathOfMatrix = "data/gene_presence_mx-big-DATorder.tsv";

			final String dirnameOfPrecalculatedSamples = args[4].trim();
			// final String dirnameOfPrecalculatedSamples =
			final String outDirName = args[5].trim();
			// final String outFileDirName = "data-processed/Maddison-resut-java";
			// "./linked_dirs/generatedSamples";

			/*
			 * blakNodeConfiguration="01" eseten a gain esemenyek adjak a fekete node-setet
			 * blakNodeConfiguration="11" eseten a present esemenyek,
			 * blakNodeConfiguration="01+11" eseten a ketto unioja, es minden ertelmes
			 * kombinacio hasznalhato
			 */
			final String blakNodeConfiguration = (args.length == 7) ? args[6].trim() : "01+11+10"; // az utolso
																									// parameternek
																									// default
																									// erteke van

			// -----------------------------

			System.out.println("Loading 0/1 matrix");
			IntMatrix x = new IntMatrix(); // ?? ebben van a gének pesence/absence állapota ??
			x.initFromTsv(new File(filepathOfMatrix));
			System.out.println("Loading 0/1 matrix- done");

			// ------------
			System.out.println("Loading tree.");
			DoubleArrayTree tree;
			if ("root".equalsIgnoreCase(rootNodeName)) {
				tree = new DoubleArrayTree(new File(treeNwkFilename));

			} else {
				tree = new DoubleArrayTree(new File(treeNwkFilename), rootNodeName);
			}
			System.out.println("Loading tree. - done");
			// ------------

			boolean ok = Arrays.equals(tree.getNames(), x.colnames);
			if (!ok) {
				throw new IllegalStateException(
						"A fa és a 0/1 tablazat nevei, vagy a nevek sorrendje nem egyezik meg.");
			}

			File precalculatedSamplesDir = new File(dirnameOfPrecalculatedSamples);

			File outDir = new File(outDirName);
			outDir.mkdirs();

			//////////////////////
			// kiszamolom elore melyek leszenk a feke regiok
			// (kb 15 sec amíg kiszamolja eg 5000 geén×365 fa-node meteru matrixra
			final boolean[][] blackMarkMatrix = calculateBlackNodeMatrix(tree, x, blakNodeConfiguration);

			int[] allTaskRows = new int[x.nrow];
			for (int i = 0; i < x.nrow; i++) {
				allTaskRows[i] = i;
			}

			// megkeverem a sorrendet
			Random random = new Random();
			for (int i = 0; i < allTaskRows.length - 1; i++) {
				int i2 = i + random.nextInt(allTaskRows.length - i);

				// megcserélem a tomb i-edik es i2-edik elemét
				int v = allTaskRows[i];
				allTaskRows[i] = allTaskRows[i2];
				allTaskRows[i2] = v;
			}

			if (N == 1) {
				System.out.println("Task order: " + Arrays.toString(allTaskRows));
				// nem kell parhuzamositani
				process(tree, x, blackMarkMatrix, allTaskRows, precalculatedSamplesDir, new File(outDir, "result.csv"));
				System.out.println("Done");
			} else {

				int n1 = allTaskRows.length / N;
				int n0 = allTaskRows.length - n1 * (N - 1);
				Thread[] workerArray = new Thread[N];

				int head = 0;
				for (int i = 0; i < N; i++) {
					int n = (i == 0 ? n0 : n1);
					int[] taskRows = new int[n];
					System.arraycopy(allTaskRows, head, taskRows, 0, n);
					head += n;

					String is = "000" + i;
					is = is.substring(is.length() - 3);
					File outFile = new File(outDir, "result-part_" + is + ".csv");

					workerArray[i] = new Thread(
							new MyWorker(tree, x, blackMarkMatrix, taskRows, precalculatedSamplesDir, outFile));
					workerArray[i].start();
					workerArray[i].setPriority(Thread.MIN_PRIORITY);
				}
				assert (head == allTaskRows.length);

				for (int i = 0; i < N; i++) {
					try {
						workerArray[i].join();
					} catch (InterruptedException e) {
						e.printStackTrace();
					}
				}
				System.out.println("All Done.");

			}
		}
	}

	private static class MyWorker implements Runnable {
		final DoubleArrayTree tree;
		final IntMatrix x;
		final int[] taskRows;
		final File precalculatedSamplesDir;
		final File outFile;
		final boolean[][] blackMarkMatrix;

		public MyWorker(DoubleArrayTree tree, IntMatrix x, boolean[][] blackMarkMatrix, int[] taskRows,
				File precalculatedSamplesDir, File outFile) {
			super();
			this.tree = tree;
			this.x = x;
			this.taskRows = taskRows;
			this.precalculatedSamplesDir = precalculatedSamplesDir;
			this.outFile = outFile;
			this.blackMarkMatrix = blackMarkMatrix;
		}

		@Override
		public void run() {
			try {
				process(tree, x, blackMarkMatrix, taskRows, precalculatedSamplesDir, outFile);
			} catch (IOException e) {
				e.printStackTrace();
			}

		}
	}

	/**
	 * Ebben a process() fuggvenben van a lenyeg
	 * 
	 * 
	 * @param tree                    az a filogenetikai fa
	 * @param x                       ebben vannak a gének jelenlétét jelző 0/1
	 *                                értékek
	 * @param blakNodeConfiguration   megadja, hogy a background gének alapján
	 *                                hogyan hatarozzuk meg mik a feket node-ok
	 *                                (vagyis elek)
	 * @param taskRows                Az x-ben benne van az osszes gen. Ez megmondja
	 *                                melyik genekre kell megcsinalni a szamitast
	 *                                (mely sorokra a matrixbol)
	 * @param precalculatedSamplesDir ebbol a konyvtarbol tolthetoek be a bootstrap
	 *                                mintak
	 * @param outFile                 output file, amibe az eredmenyt irjuk.
	 *                                szalankent kulon output fileok vanna.
	 * 
	 * @throws FileNotFoundException
	 * @throws IOException
	 */

	private static void process(DoubleArrayTree tree, IntMatrix x, boolean[][] blackMarkMatrix, int[] taskRows,
			File precalculatedSamplesDir, final File outFile) throws FileNotFoundException, IOException {

		// fekete él az, ahol a G2~PRESENCE
		// G1 GAIN és LOSS számai alapkján generál bootsrap mintákat a fán lehetséges
		// gén 0/1-re

		BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outFile)));

		try {

			writer.write(HEADER);
			for (int iTask = 0; iTask < taskRows.length; iTask++) {

				int rowI1 = taskRows[iTask];

				String cluster1 = x.rownames[rowI1]; // ez a gén a TESTED
				System.out.println("TESTED GENE:" + cluster1);// DEBUG
				System.out.println("STARTING ROW [" + rowI1 + "/" + x.nrow + "  -- " + cluster1 + "  ~  " + iTask + "/"
						+ taskRows.length + "]");

				int[] gene1Values = x.getRow(rowI1);
				int[] cntArrG1 = DATUtils.count_loss_and_gain_eventses(gene1Values, tree, tree.getRoot());

				Precalculated precalculated = new MaddisonSimulationTestDATFromPreparedData.Precalculated(
						cntArrG1[GAIN], cntArrG1[LOSS], tree, precalculatedSamplesDir);
				int bootstrapSampleSize = precalculated.labelingList1.size();

				for (int rowI2 = 0; rowI2 < x.nrow; rowI2++) {
					// 	ez a for megy végig a BACGROUND geneken

					String cluster2 = x.rownames[rowI2]; // Ez a gene BACKGROUND
					// System.out.println(cluster1 + "-" + cluster2); //DEBUG
					int[] gene2Values = x.getRow(rowI2);

					boolean[] black_idx_based_on_g2 = blackMarkMatrix[rowI2]; // ez mar elore ki van szamitva


					int blackCnt = sum(black_idx_based_on_g2);
					double blackRatio = ((double) blackCnt) / black_idx_based_on_g2.length;

//					System.out.println("blackRatio="+blackRatio); //DEBUG
					int[] actualCnts = DATUtils.count_loss_and_gain_eventses_in_black_area(gene1Values, tree,
							black_idx_based_on_g2);

					// contingency matrix represented ina vector
					int[] contingency_vector1 = DATUtils.count_loss_and_gain_eventses(gene1Values, gene2Values, tree,
							tree.getRoot());

					String tableAsString;
					String resultPartStringOfPValues;
					if (bootstrapSampleSize > 0) {
						int[] table = precalculated.calulateTable(black_idx_based_on_g2);

						assert (cntArrG1[GAIN] + 1) * (cntArrG1[LOSS] + 1) == table.length;

						StringBuilder sb = new StringBuilder();
						for (int i = 0; i < table.length; i++) {
							sb.append(table[i]).append('_');
						}
						sb.setLength(sb.length() - 1);
						tableAsString = sb.toString();

						assert (sum(table) == bootstrapSampleSize);

						double[] probablibilityTable = divide(table, bootstrapSampleSize);
						double[] loss_probablibility_vector = rowSum(probablibilityTable, /* nrow */cntArrG1[LOSS] + 1,
								/* ncol */cntArrG1[GAIN] + 1);
						double[] gain_probablibility_vector = colSum(probablibilityTable, /* nrow */cntArrG1[LOSS] + 1,
								/* ncol */cntArrG1[GAIN] + 1);
						assert (loss_probablibility_vector.length == cntArrG1[LOSS] + 1);
						assert (gain_probablibility_vector.length == cntArrG1[GAIN] + 1);
						assert (equalsDouble(sum(probablibilityTable), 1.0));
						assert (equalsDouble(sum(loss_probablibility_vector), 1.0));
						assert (equalsDouble(sum(gain_probablibility_vector), 1.0));

						double expedted_number_of_black_gains = expectedValue(gain_probablibility_vector);
						double expedted_number_of_black_loss = expectedValue(loss_probablibility_vector);
						double estimated_median_of_black_gains = my_median(gain_probablibility_vector);
						double estimated_median_of_black_loss = my_median(loss_probablibility_vector);

						double p_value_less_gains_then_expected = partialSum(gain_probablibility_vector, 0,
								actualCnts[GAIN]);
						double p_value_more_gains_then_expected = partialSum(gain_probablibility_vector,
								actualCnts[GAIN], cntArrG1[GAIN]);
						double p_value_less_loss_then_expected = partialSum(loss_probablibility_vector, 0,
								actualCnts[LOSS]);
						double p_value_more_loss_then_expected = partialSum(loss_probablibility_vector,
								actualCnts[LOSS], cntArrG1[LOSS]);

						resultPartStringOfPValues = "" + expedted_number_of_black_gains + ","
								+ expedted_number_of_black_loss + "," + estimated_median_of_black_gains + ","
								+ estimated_median_of_black_loss + "," + p_value_less_gains_then_expected + ","
								+ p_value_more_gains_then_expected + "," + p_value_less_loss_then_expected + ","
								+ p_value_more_loss_then_expected + ",";

					} else {
						tableAsString = "";
						resultPartStringOfPValues = "NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,";
					}

					writer.write(cluster1);
					writer.write(',');
					writer.write(cluster2);
					writer.write(',');

					writer.write(String.valueOf(blackRatio));
					writer.write(',');

					writer.write(String.valueOf(bootstrapSampleSize));
					writer.write(',');

					//

					writeAllComaSeparated(writer, cntArrG1); // G1
					writeAllComaSeparated(writer, actualCnts); // black
					writeAllComaSeparated(writer, contingency_vector1); // contingencia tabla vektor formában
					writer.write(resultPartStringOfPValues);
					writer.write(tableAsString);

					writer.write('\n');
				}
			}
		} finally {
			writer.close();
		}
	}

	private static boolean[][] calculateBlackNodeMatrix(DoubleArrayTree tree, IntMatrix x,
			String blakNodeConfiguration) {

		boolean[][] blackMarkMatrix = new boolean[x.nrow][];
		for (int rowI2 = 0; rowI2 < x.nrow; rowI2++) {

			// String cluster2 = x.rownames[rowI2]; // Ez a gene BACKGROUND
			// System.out.println(cluster1 + "-" + cluster2); //DEBUG
			int[] gene2Values = x.getRow(rowI2);

			boolean[] black_idx_based_on_g2 = calculateBlackNodes(tree, blakNodeConfiguration, gene2Values);
			blackMarkMatrix[rowI2] = black_idx_based_on_g2;
		}

		return blackMarkMatrix;
	}

	private static boolean[] calculateBlackNodes(DoubleArrayTree tree, String blakNodeConfiguration,
			int[] gene2Values) {
		// -------------------------
		// Osszeallitja azon node-ok indexeit, amiket feketere kell szinezni.
		// A blakNodeConfiguration string alapjan donti el,
		// hogy a G2 gén GAIN/PRESENCE/LOSS/ABSENCE esemeyeibol melyek legyenek feketek
		boolean[] black_idx_based_on_g2 = null;

//					System.out.println("BACKGROUND STATE cnts"+ Arrays.toString(DATUtils.count_gene_presence_absence(gene2Values, tree,tree.getRoot())));// DEBUG
//					System.out.println("BACKGROUND ACTION cnts"+ Arrays.toString(DATUtils.count_loss_and_gain_eventses(gene2Values, tree,tree.getRoot())));// DEBUG 
		if (blakNodeConfiguration.contains("01")) {
			boolean[] black_idx_tmp = DATUtils.selectNodes(gene2Values, tree, 0, 1);
			black_idx_based_on_g2 = (black_idx_based_on_g2 == null) ? black_idx_tmp
					: DATUtils.orOperatorOfBooleanArrays(black_idx_based_on_g2, black_idx_tmp);
		}
		if (blakNodeConfiguration.contains("10")) {
			boolean[] black_idx_tmp = DATUtils.selectNodes(gene2Values, tree, 1, 0);
			black_idx_based_on_g2 = (black_idx_based_on_g2 == null) ? black_idx_tmp
					: DATUtils.orOperatorOfBooleanArrays(black_idx_based_on_g2, black_idx_tmp);
		}
		if (blakNodeConfiguration.contains("00")) {
			boolean[] black_idx_tmp = DATUtils.selectNodes(gene2Values, tree, 0, 0);
			black_idx_based_on_g2 = (black_idx_based_on_g2 == null) ? black_idx_tmp
					: DATUtils.orOperatorOfBooleanArrays(black_idx_based_on_g2, black_idx_tmp);
		}
		if (blakNodeConfiguration.contains("11")) {
			boolean[] black_idx_tmp = DATUtils.selectNodes(gene2Values, tree, 1, 1);
			black_idx_based_on_g2 = (black_idx_based_on_g2 == null) ? black_idx_tmp
					: DATUtils.orOperatorOfBooleanArrays(black_idx_based_on_g2, black_idx_tmp);
		}

		if (black_idx_based_on_g2 == null) {
			throw new IllegalArgumentException("String blakNodeConfiguration==\"" + blakNodeConfiguration + "\"");
		}
		return black_idx_based_on_g2;
	}

	private static final void writeAllComaSeparated(Writer writer, int[] arr) throws IOException {
		for (int i = 0; i < arr.length; i++) {
			writer.write(String.valueOf(arr[i]));
			writer.write(',');
		}

	}

	private static final int sum(boolean[] arr) {
		int sum = 0;
		for (int i = 0; i < arr.length; i++) {
			if (arr[i]) {
				sum++;
			}
		}
		return sum;
	}

	/** pontonkenti osztas konstanssal ugy, hogy int[]-bol double[]-ot csinal */
	private static final double[] divide(int[] arr, int x) {
		double xd = x;
		double[] res = new double[arr.length];
		for (int i = 0; i < arr.length; i++) {
			res[i] = ((double) arr[i]) / xd;
		}
		return res;

	}

	private static final double sum(double[] arr) {
		double sum = 0;
		for (int i = 0; i < arr.length; i++) {
			sum += arr[i];
		}
		return sum;
	}

	private static final int sum(int[] arr) {
		int sum = 0;
		for (int i = 0; i < arr.length; i++) {
			sum += arr[i];
		}
		return sum;
	}

	/**
	 * összegzi az arr[] tömb egy részét, a from és to indexet iinclusive-ankezeli
	 */
	private static final double partialSum(double[] arr, int from, int to) {
		assert (from <= to);
		double sum = 0;
		for (int i = from; i <= to; i++) {
			sum += arr[i];
		}
		return sum;
	}

	/**
	 * Ugy tekintem mintha az arr-ban levo vektor egy matrix lenne az ncol×nrov
	 * méretekkel, és a sorainak az összget adom vissza egy double[nrow] meretu
	 * vektiorban
	 * 
	 * A tablazatot ugy kepzelem, hogy oszlop-folytonosan van ábrázolva, tehát az
	 * első néhány elem az arr-ban az első oszlopban van. Ha a kovetkezo sorra
	 * akarunk lépni, akkor nrow-val kel növelni az indexet az arr[]-ban
	 */
	private static final double[] rowSum(double[] arr, int nrow, int ncol) {
		assert (nrow * ncol == arr.length);

		double[] res = new double[nrow];
		for (int i = 0; i < ncol; i++) // balrol-jobbra
		{
			int i_X_nrow = i * nrow;
			for (int j = 0; j < nrow; j++) // fentrol-le
			{
				res[j] += arr[i_X_nrow + j];
			}

		}
		return res;

	}

	private static final double[] colSum(double[] arr, int nrow, int ncol) {
		assert (nrow * ncol == arr.length);

		double[] res = new double[ncol];
		for (int i = 0; i < ncol; i++) {
			int i_X_nrow = i * nrow;
			for (int j = 0; j < nrow; j++) {
				res[i] += arr[i_X_nrow + j];
			}

		}
		return res;

	}

	private static final double my_median(double[] prob_vektor) {
		assert (Math.abs(sum(prob_vektor) - 1) < 0.0001); // sum(prob_vektor) == 1 , csak a szamitasi pontossag miatt
		// bonyolult

		double p = 0;
		double sum = 0;
		for (int i = 0; i < prob_vektor.length; i++) {
			sum += prob_vektor[i];
			if (sum > 0.5 & p < 0.5) {
				return i;
			}
			if (sum == 0.5) {
				return i + 0.5;
			}
			p = sum;
		}
		throw new RuntimeException("median szamitasi hiba+ " + Arrays.toString(prob_vektor));

		// # prob_vektor<-rep(0.25,4)
		// prob_vektor<-cumsum(prob_vektor)
		// prob_vektor_lag<-lag(prob_vektor)
		// prob_vektor_lag[1]<-0
		// idx<-prob_vektor_lag<=0.5 & prob_vektor>=0.5
		// values<-which(idx)-1
		// return( mean(values))
	}

	private static final double expectedValue(double[] prob_vektor) {
		double r = 0;
		for (int i = 0; i < prob_vektor.length; i++) {
			r += i * prob_vektor[i];
		}
		return r;
	}

	private static final boolean equalsDouble(double x, double y) {
		return Math.abs(x - y) < 0.0001;
	}

	/** test the rowSum anc ColSum functioh */
	public static void main3(String[] args) {
		double[] x = new double[] { 1, 2, 3, 4, 5, 6 };
		double[] y = rowSum(x, 2, 3); // rowSum(double[] arr, int nrow, int ncol)

		System.out.println(Arrays.toString(y));

		y = colSum(x, 2, 3);
		System.out.println(Arrays.toString(y));

		x = new double[] { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };

		y = rowSum(x, 4, 3);
		System.out.println(Arrays.toString(y));

		y = colSum(x, 4, 3);
		System.out.println(Arrays.toString(y));

	}

}
