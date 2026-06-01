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



/**
 * Run the simulation-based Maddison test for all foreground/background gene pairs.
 *
 * <p>
 * This class performs the main Maddison-test calculation in the pipeline. For
 * each foreground gene and each background gene, it tests whether gain or loss
 * events of the foreground gene are enriched in the black region of the
 * phylogenetic tree defined by the background gene.
 * </p>
 *
 * <p>
 * The input gene presence/absence matrix is expected to be binary and ordered
 * according to the node order of the Java {@code DoubleArrayTree}
 * representation. Rows correspond to genes or orthogroups, and columns
 * correspond to tree nodes:
 * </p>
 *
 * <pre>
 * 0 = gene absent
 * 1 = gene present
 * </pre>
 *
 * <p>
 * Along each parent-to-child tree edge, the program distinguishes four
 * transition categories:
 * </p>
 *
 * <pre>
 * 0 -> 0  absent
 * 0 -> 1  gain
 * 1 -> 0  loss
 * 1 -> 1  present
 * </pre>
 *
 * <p>
 * In the output table, the foreground gene is stored in the {@code cluster1}
 * column and the background gene is stored in the {@code cluster2} column.
 * These column names are kept for compatibility with downstream analysis steps.
 * </p>
 *
 * <h2>Black-region definition</h2>
 *
 * <p>
 * The black region is defined separately for each background gene. Events are
 * biologically associated with tree edges, but in the implementation each edge
 * is represented by the index of its child node. Therefore, a black region is
 * stored as a {@code boolean[]} mask over tree node indices.
 * </p>
 *
 * <p>
 * The optional {@code blackNodeConfiguration} command-line argument specifies
 * which background-gene parent-to-child transitions are included in the black
 * region:
 * </p>
 *
 * <pre>
 * 00  background gene absent at both endpoints
 * 01  background-gene gain
 * 10  background-gene loss
 * 11  background gene present at both endpoints
 * </pre>
 *
 * <p>
 * Multiple categories can be combined with {@code +}, for example
 * {@code 01+11+10}. In the current biological pipeline, the main use case is
 * {@code 11}, meaning that an edge is black when the background gene is present
 * at both endpoints.
 * </p>
 *
 * <h2>Simulation model</h2>
 *
 * <p>
 * For a foreground gene, let {@code G} be the total number of gains and
 * {@code L} the total number of losses on the tree. The program loads
 * pregenerated random tree labelings with the same total {@code (G, L)} counts.
 * These labelings are produced by {@code PrepareMaddisonSimulation}.
 * </p>
 *
 * <p>
 * For each background gene, the program counts in each simulated labeling:
 * </p>
 *
 * <pre>
 * x = number of gains in the black region
 * y = number of losses in the black region
 * </pre>
 *
 * <p>
 * This gives a Monte Carlo estimate of the two-dimensional discrete
 * distribution of black-region gain and loss counts conditional on the total
 * foreground-gene gain/loss counts and on the background-gene black region.
 * The estimated frequency table is written to the
 * {@code table_underscore_separated} output column.
 * </p>
 *
 * <p>
 * From this two-dimensional frequency table, the program also derives
 * one-dimensional gain and loss distributions. These are used to estimate the
 * expected number and median number of black-region gains and losses, and to
 * calculate one-sided p-values:
 * </p>
 *
 * <pre>
 * p_value_less_gains_than_expected
 * p_value_more_gains_than_expected
 * p_value_less_loss_than_expected
 * p_value_more_loss_than_expected
 * </pre>
 *
 * <h2>Additional output</h2>
 *
 * <p>
 * The output is a CSV table. Each row corresponds to one
 * {@code (foreground gene, background gene)} pair. Besides the Maddison p-values,
 * the table also contains observed transition counts, simulation summary
 * statistics, and a 4x4 contingency matrix describing the joint transition
 * categories of the foreground and background genes on tree edges.
 * </p>
 *
 * <p>
 * The contingency-matrix columns use the {@code CM_} prefix. The four transition
 * categories of the foreground gene are crossed with the four transition
 * categories of the background gene, giving sixteen columns such as:
 * </p>
 *
 * <pre>
 * CM_G1_absent_G2_absent
 * CM_G1_gain_G2_absent
 * CM_G1_loss_G2_present
 * CM_G1_present_G2_present
 * </pre>
 *
 * <h2>Parallel execution</h2>
 *
 * <p>
 * Parallelization is implemented by splitting foreground genes among worker
 * threads. Each worker processes its assigned foreground genes against all
 * background genes and writes a separate CSV result part. A later R pipeline
 * step concatenates these result files.
 * </p>
 *
 * <h2>Missing simulation samples</h2>
 *
 * <p>
 * If no pregenerated simulation samples are available for a foreground gene's
 * total {@code (G, L)} count, the program cannot estimate the null distribution.
 * In that case, {@code sample_size_of_simulation} is set to {@code 0}, the
 * frequency-table column is empty, and p-values and summary statistics are
 * written as {@code NaN}.
 * </p>
 */
public class DoMaddisonTestForAllPairs {
	static final int ABSENT = 0;
	static final int GAIN = 1;
	static final int LOSS = 2;
	static final int PRESENT = 3;

	private static final String HEADER_PART_OF_P_VALUES = "expected_number_of_black_gains,expected_number_of_black_loss,estimated_median_of_black_gains,estimated_median_of_black_loss,"
			+ "p_value_less_gains_than_expected,p_value_more_gains_than_expected,p_value_less_loss_than_expected,p_value_more_loss_than_expected,";
	private static final String HEADER = "cluster1,cluster2,black_area_ratio,sample_size_of_simulation," + //
			"G1_absent,G1_gain,G1_loss,G1_present," + //
			"actual_black_absent,actual_black_gain,actual_black_loss,actual_black_present," + //
			"CM_G1_absent_G2_absent,CM_G1_gain_G2_absent,CM_G1_loss_G2_absent,CM_G1_present_G2_absent,CM_G1_absent_G2_gain,CM_G1_gain_G2_gain,CM_G1_loss_G2_gain,CM_G1_present_G2_gain,CM_G1_absent_G2_loss,CM_G1_gain_G2_loss,CM_G1_loss_G2_loss,CM_G1_present_G2_loss,CM_G1_absent_G2_present,CM_G1_gain_G2_present,CM_G1_loss_G2_present,CM_G1_present_G2_present,"
			+ //
			HEADER_PART_OF_P_VALUES + //
			"table_underscore_separated" + //
			"\n";

	public static void main(String[] args) throws IOException {

		System.out.println(Arrays.toString(args));
		if (args.length != 6 & args.length != 7) {
			  System.err.println("Expected 6 or 7 command-line arguments:\n"
			            + "  [1] numOfWorkers\n"
			            + "      Number of parallel worker threads.\n"
			            + "  [2] treeNwkFilename\n"
			            + "      Path to the input phylogenetic tree in Newick format.\n"
			            + "  [3] rootNodeName\n"
			            + "      Name of the root node of the subtree to analyse. Use 'root' for the full tree.\n"
			            + "  [4] genePresenceMatrixFilename\n"
			            + "      Path to the DAT-ordered binary gene presence/absence matrix.\n"
			            + "  [5] dirnameOfPrecalculatedSamples\n"
			            + "      Directory containing pregenerated random tree labelings.\n"
			            + "      (This directory is the result of PrepareMaddisonSimulation.java .)\n"
			            + "  [6] outDirName\n"
			            + "      Output directory where result CSV files will be written.\n"
			            + "  [7] blackNodeConfiguration\n"
			            + "      Optional definition of the black region based on background-gene transitions,\n"
			            + "      for example '11'. Default: '01+11+10'.\n");
			    System.exit(1);
		} else {
			final int N = Integer.parseInt(args[0]); // Number of parallel worker threads.
			String treeNwkFilename = args[1].trim();
			String rootNodeName = args[2].trim();
			final String filepathOfMatrix = args[3].trim();
			final String dirnameOfPrecalculatedSamples = args[4].trim();
			final String outDirName = args[5].trim();

			/*
			 * The blackNodeConfiguration parameter defines which background-gene
			 * parent-to-child transitions are included in the black region.
			 *
			 * Examples:
			 *   "01"       include background-gene gain edges
			 *   "11"       include edges where the background gene is present at both endpoints
			 *   "01+11"    include the union of both categories
			 *
			 * Any meaningful combination of "00", "01", "10", and "11" can be used.
			 */
			final String blackNodeConfiguration = (args.length == 7) ? args[6].trim() : "01+11+10"; // default value
			// -----------------------------

			System.out.println("Loading 0/1 matrix");
			IntMatrix x = new IntMatrix();  // Binary gene presence/absence matrix.
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
				        "Tree node names and matrix column names do not match, or they are not in the same order.");
			}

			File precalculatedSamplesDir = new File(dirnameOfPrecalculatedSamples);

			File outDir = new File(outDirName);
			outDir.mkdirs();

			//////////////////////
			// Precompute the black-region mask for each background gene.
			// This avoids recalculating the black region for every foreground/background pair.
			final boolean[][] blackMarkMatrix = calculateBlackNodeMatrix(tree, x, blackNodeConfiguration);

			int[] allTaskRows = new int[x.nrow];
			for (int i = 0; i < x.nrow; i++) {
				allTaskRows[i] = i;
			}

			// Shuffle the foreground-gene processing order.
			Random random = new Random();
			for (int i = 0; i < allTaskRows.length - 1; i++) {
				int i2 = i + random.nextInt(allTaskRows.length - i);

				// Swap elements i and i2.
				int v = allTaskRows[i];
				allTaskRows[i] = allTaskRows[i2];
				allTaskRows[i2] = v;
			}

			if (N == 1) {
				// No parallelization is needed when only one worker is used.
				System.out.println("Task order: " + Arrays.toString(allTaskRows));
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
	 * Run the Maddison simulation test for the selected foreground genes.
	 *
	 * Each row index in {@code taskRows} identifies one foreground gene in the
	 * binary gene presence/absence matrix. For each foreground gene, the method
	 * iterates over all possible background genes. The background gene defines the
	 * black region of the tree, and the foreground gene provides the gain/loss
	 * events to be tested for enrichment in that region.
	 *
	 * Each worker thread calls this method with a different subset of foreground
	 * genes and writes its results to a separate CSV file.
	 *
	 * @param tree
	 *            phylogenetic tree represented as a DoubleArrayTree
	 * @param x
	 *            DAT-ordered binary gene presence/absence matrix
	 * @param blackMarkMatrix
	 *            precomputed black-region masks, one per background gene
	 * @param taskRows
	 *            row indices of the foreground genes assigned to this worker (x contains all the genes, but this worker has to do the indexed ones)
	 * @param precalculatedSamplesDir
	 *            directory containing pregenerated random tree labelings
	 * @param outFile
	 *            output CSV file written by this worker
	 *
	 * @throws FileNotFoundException
	 *             if an input or output file cannot be opened
	 * @throws IOException
	 *             if reading or writing fails
	 */
	private static void process(DoubleArrayTree tree, IntMatrix x, boolean[][] blackMarkMatrix, int[] taskRows,
			File precalculatedSamplesDir, final File outFile) throws FileNotFoundException, IOException {

		// The black region is defined by the background gene.
		// The foreground gene provides the observed gain/loss counts and determines
		// which pregenerated random labelings are loaded.
		
		BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outFile)));

		try {

			writer.write(HEADER);
			for (int iTask = 0; iTask < taskRows.length; iTask++) {

				int rowI1 = taskRows[iTask];

				String cluster1 = x.rownames[rowI1]; // Foreground/tested gene.
				System.out.println("TESTED GENE:" + cluster1);// DEBUG
				System.out.println("STARTING ROW [" + rowI1 + "/" + x.nrow + "  -- " + cluster1 + "  ~  " + iTask + "/"
						+ taskRows.length + "]");

				int[] gene1Values = x.getRow(rowI1);
				int[] cntArrG1 = DATUtils.count_loss_and_gain_eventses(gene1Values, tree, tree.getRoot());

				Precalculated precalculated = new MaddisonSimulationTestDATFromPreparedData.Precalculated(
						cntArrG1[GAIN], cntArrG1[LOSS], tree, precalculatedSamplesDir);
				int bootstrapSampleSize = precalculated.labelingList1.size();

				for (int rowI2 = 0; rowI2 < x.nrow; rowI2++) {
					// Iterate over all background genes.

					String cluster2 = x.rownames[rowI2]; // Background gene.
					int[] gene2Values = x.getRow(rowI2);

					boolean[] black_idx_based_on_g2 = blackMarkMatrix[rowI2]; // Precomputed black-region mask for this background gene.


					int blackCnt = sum(black_idx_based_on_g2);
					double blackRatio = ((double) blackCnt) / black_idx_based_on_g2.length;

					//	System.out.println("blackRatio="+blackRatio); //DEBUG
					int[] actualCnts = DATUtils.count_loss_and_gain_eventses_in_black_area(gene1Values, tree,
							black_idx_based_on_g2);

					// Contingency matrix represented as a vector.
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

						double expected_number_of_black_gains = expectedValue(gain_probablibility_vector);
						double expected_number_of_black_loss = expectedValue(loss_probablibility_vector);
						double estimated_median_of_black_gains = my_median(gain_probablibility_vector);
						double estimated_median_of_black_loss = my_median(loss_probablibility_vector);

						double p_value_less_gains_than_expected = partialSum(gain_probablibility_vector, 0,
								actualCnts[GAIN]);
						double p_value_more_gains_than_expected = partialSum(gain_probablibility_vector,
								actualCnts[GAIN], cntArrG1[GAIN]);
						double p_value_less_loss_than_expected = partialSum(loss_probablibility_vector, 0,
								actualCnts[LOSS]);
						double p_value_more_loss_than_expected = partialSum(loss_probablibility_vector,
								actualCnts[LOSS], cntArrG1[LOSS]);

						resultPartStringOfPValues = "" + expected_number_of_black_gains + ","
								+ expected_number_of_black_loss + "," + estimated_median_of_black_gains + ","
								+ estimated_median_of_black_loss + "," + p_value_less_gains_than_expected + ","
								+ p_value_more_gains_than_expected + "," + p_value_less_loss_than_expected + ","
								+ p_value_more_loss_than_expected + ",";

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
					writeAllComaSeparated(writer, contingency_vector1); // Contingency table in vector form.
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
			String blackNodeConfiguration) {

		boolean[][] blackMarkMatrix = new boolean[x.nrow][];
		for (int rowI2 = 0; rowI2 < x.nrow; rowI2++) {

			// String cluster2 = x.rownames[rowI2]; // Background gene.
			// System.out.println(cluster1 + "-" + cluster2); //DEBUG
			int[] gene2Values = x.getRow(rowI2);

			boolean[] black_idx_based_on_g2 = calculateBlackNodes(tree, blackNodeConfiguration, gene2Values);
			blackMarkMatrix[rowI2] = black_idx_based_on_g2;
		}

		return blackMarkMatrix;
	}

	private static boolean[] calculateBlackNodes(DoubleArrayTree tree, String blackNodeConfiguration,
			int[] gene2Values) {
		// -------------------------
		// Build the boolean node mask that marks the black region.
		// The blackNodeConfiguration string determines which background-gene (G2)
		// parent-to-child transition categories (GAIN/PRESENCE/LOSS/ABSENCE) are included in the black region.
		boolean[] black_idx_based_on_g2 = null;

		//	System.out.println("BACKGROUND STATE cnts"+ Arrays.toString(DATUtils.count_gene_presence_absence(gene2Values, tree,tree.getRoot())));// DEBUG
		//	System.out.println("BACKGROUND ACTION cnts"+ Arrays.toString(DATUtils.count_loss_and_gain_eventses(gene2Values, tree,tree.getRoot())));// DEBUG 
		if (blackNodeConfiguration.contains("01")) {
			boolean[] black_idx_tmp = DATUtils.selectNodes(gene2Values, tree, 0, 1);
			black_idx_based_on_g2 = (black_idx_based_on_g2 == null) ? black_idx_tmp
					: DATUtils.orOperatorOfBooleanArrays(black_idx_based_on_g2, black_idx_tmp);
		}
		if (blackNodeConfiguration.contains("10")) {
			boolean[] black_idx_tmp = DATUtils.selectNodes(gene2Values, tree, 1, 0);
			black_idx_based_on_g2 = (black_idx_based_on_g2 == null) ? black_idx_tmp
					: DATUtils.orOperatorOfBooleanArrays(black_idx_based_on_g2, black_idx_tmp);
		}
		if (blackNodeConfiguration.contains("00")) {
			boolean[] black_idx_tmp = DATUtils.selectNodes(gene2Values, tree, 0, 0);
			black_idx_based_on_g2 = (black_idx_based_on_g2 == null) ? black_idx_tmp
					: DATUtils.orOperatorOfBooleanArrays(black_idx_based_on_g2, black_idx_tmp);
		}
		if (blackNodeConfiguration.contains("11")) {
			boolean[] black_idx_tmp = DATUtils.selectNodes(gene2Values, tree, 1, 1);
			black_idx_based_on_g2 = (black_idx_based_on_g2 == null) ? black_idx_tmp
					: DATUtils.orOperatorOfBooleanArrays(black_idx_based_on_g2, black_idx_tmp);
		}

		if (black_idx_based_on_g2 == null) {
			throw new IllegalArgumentException("Invalid blackNodeConfiguration: \"" + blackNodeConfiguration + "\"");
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

	
	/** Divide each element by a constant and convert the int array to a double array. */
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
	 * Sum a range of an array, treating both {@code from} and {@code to} as inclusive indices.
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
	 * Treat a vector as a column-major matrix and return its row sums.
	 *
	 * The vector is interpreted as a matrix with {@code nrow} rows and {@code ncol}
	 * columns. Consecutive elements in the vector belong to the same column. To
	 * move to the next column, the index is increased by {@code nrow}.
	 */
	private static final double[] rowSum(double[] arr, int nrow, int ncol) {
		assert (nrow * ncol == arr.length);

		double[] res = new double[nrow];
		for (int i = 0; i < ncol; i++) // left to right
		{
			int i_X_nrow = i * nrow;
			for (int j = 0; j < nrow; j++) // top to bottom
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
		assert (Math.abs(sum(prob_vektor) - 1) < 0.0001); // sum(prob_vektor) == 1 , and  allow a small numerical tolerance.

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
		throw new RuntimeException("Median calculation failed: " + Arrays.toString(prob_vektor));
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
	
}
