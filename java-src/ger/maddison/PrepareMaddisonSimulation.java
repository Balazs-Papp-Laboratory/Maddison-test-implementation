package ger.maddison;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Random;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.ZipEntry;
import java.util.zip.ZipInputStream;
import java.util.zip.ZipOutputStream;

import ger.maddison.MaddisonSimulationTestDAT.LabelingGandAndLoss;
import ger.maddison.MyTaskList.Row;
import ger.tree.DoubleArrayTree;
import ger.tree.IntegerStack;


/**
 * Pregenerate random tree labelings for the simulation-based Maddison test.
 *
 * This program works together with
 * {@code script-03-create_simulation_task_table.R}. The R script reads the
 * per-gene gain/loss counts produced by the previous pipeline step and creates
 * a task table. Each row of that table describes one group of simulation
 * samples to generate: a fixed total number of changes
 * ({@code gain_plus_loss}) and the set of gain counts that are needed for that
 * total.
 *
 * The Java program reads this task table and generates random binary
 * presence/absence labelings on a fixed phylogenetic tree. A labeling assigns
 * one binary state to each tree node:
 *
 * <pre>
 * 0 = gene absent
 * 1 = gene present
 * </pre>
 *
 * For each task, the program fixes the total number of state changes on the
 * tree. It then repeatedly selects random change positions, propagates binary
 * states along the tree, and counts how many of the induced changes are gains
 * and how many are losses:
 *
 * <pre>
 * 0 -> 1  gain
 * 1 -> 0  loss
 * </pre>
 *
 * The program uses an accept/reject strategy. It can directly control the total
 * number of changes, but not the exact split between gains and losses. Therefore
 * it keeps only those random labelings whose observed gain count matches one of
 * the gain counts requested in the task table. Other generated labelings are
 * discarded.
 *
 * Grouping tasks by {@code gain_plus_loss} improves efficiency. Several
 * different {@code (gain, loss)} pairs can share the same total number of
 * changes. A single block of generated labelings can therefore contribute
 * samples to multiple requested gain/loss categories.
 *
 * Generated labelings are written to ZIP files in the output directory. The
 * program also maintains a broker file named {@code _BROKER_FILE.txt}. This file
 * records which simulation samples have already been generated and where they
 * are stored. The broker file allows long runs to be interrupted and resumed:
 * when the program is restarted, it reads the broker state and continues from
 * the previously generated samples.
 *
 * The pregenerated samples produced by this program are used by the next
 * pipeline step, which performs the actual simulation-based Maddison test and
 * calculates p-values for gene pairs.
 */
public class PrepareMaddisonSimulation {
	private static class BrokerTaskUnit {
		public int gainPlusLoss;
		// public final int n;
		public boolean[] interestingGainIdx;
		public boolean[] originalInterestingGainIdx;
		public int[] cntArray;
		public int cntTrial;

		public BrokerTaskUnit(MyTaskList.Row r) {
			super();
			this.gainPlusLoss = r.gainPlusLoss;
			this.interestingGainIdx = r.gains;
			this.originalInterestingGainIdx = Arrays.copyOf(r.gains, r.gains.length);
			this.cntArray = new int[interestingGainIdx.length];
			this.cntTrial = 0;

		}

		void addHistogramToBrokerTaskUnit(int[] histogram, int iterations) {
			final int gainPlusLoss = histogram.length - 1;
			for (int i = 0; i < histogram.length; i++) {
				if (interestingGainIdx[i]) {
					cntArray[i] += histogram[i];
					if (i != gainPlusLoss - i) {
						// Counts are mirrored because gain and loss labels are swapped when the
						// root state is changed from 1 to 0.
						cntArray[i] += histogram[gainPlusLoss - i];
					}
				}
			}

			cntTrial += iterations;
		}

		int getMinimalInterestingCount() {
			int min = Integer.MAX_VALUE;
			for (int i = 0; i < interestingGainIdx.length; i++) {
				if (interestingGainIdx[i]) {
					int x = cntArray[i];
					if (x < min) {
						min = x;
					}
				}
			}
			return min;
		}

		void updateInterestingGainIndex(int maxIteration) {
			for (int i = 0; i < cntArray.length; i++) {
				if (interestingGainIdx[i]) {
					interestingGainIdx[i] = cntArray[i] < maxIteration;
				}
			}
		}

		@Override
		public String toString() {
			StringBuilder sb = new StringBuilder();
			sb.append("(gain+loss=").append(gainPlusLoss);
			for (int i = 0; i < interestingGainIdx.length; i++) {
				if (interestingGainIdx[i]) {
					sb.append(" ");
					sb.append(i).append('[').append(cntArray[i]).append(']');
				} else if (cntArray[i] > 0) {
					sb.append(" *").append(i).append('[').append(cntArray[i]).append(']');
				}
			}
			sb.append(")\n");
			return sb.toString();
		}

		public String toString(int minIteration, int maxIteration, int maxTrial) {
			StringBuilder sb = new StringBuilder();
			sb.append("(gain+loss=").append(gainPlusLoss);
			sb.append(" cntTrial=").append(cntTrial);
			if (cntTrial >= maxTrial) {
				sb.append('#');
			}
			for (int i = 0; i < originalInterestingGainIdx.length; i++) {
				if (originalInterestingGainIdx[i]) {
					int n = cntArray[i];
					char decorator = (n >= maxIteration) ? '*' : (n >= minIteration) ? '+' : '.';
					sb.append(' ').append(i).append('[').append(n).append(decorator).append(']');
				}
			}
			sb.append(")");
			return sb.toString();
		}

	}
	/// ------------------------------------------------

	private static class TaskBroker {
		private final int minIteration; // Minimum number of accepted samples required for each target gain count.
		private final int maxIteration; // Do not store additional samples for a gain count after this many hits.
		private final int maxTrial; // Stop trying after this many generated labelings for a task.

		private final String logFileName;

		private final LinkedHashMap<Integer, BrokerTaskUnit> taskMapFree;
		private final HashMap<Integer, BrokerTaskUnit> taskMapUnderWork;
		private final ArrayList<BrokerTaskUnit> taskListFinished;

		private final Writer brokerFileWriter;

		public TaskBroker(ArrayList<Row> taskList, int minIteration, int maxIteration, int maxTrial, String logFileName,
				File brokerFile) throws IOException {
			this.minIteration = minIteration;
			this.maxIteration = maxIteration;
			this.maxTrial = maxTrial;
			this.logFileName = logFileName;
			taskMapFree = new LinkedHashMap<Integer, BrokerTaskUnit>();
			taskMapUnderWork = new HashMap<Integer, BrokerTaskUnit>();
			taskListFinished = new ArrayList<BrokerTaskUnit>();
			for (Row r : taskList) {
				taskMapFree.put(r.gainPlusLoss, new BrokerTaskUnit(r));

			}

			if (brokerFile.exists()) {
				File previousBrokerFile = new File(brokerFile.getAbsoluteFile() + ".OLD");
				brokerFile.renameTo(previousBrokerFile);

				brokerFileWriter = new OutputStreamWriter(new FileOutputStream(brokerFile), "UTF-8");
				brokerFileWriter.write("gainPlusLoss;iterations;filename;histogramArray\n");

				final Pattern lineRegexpPattern = Pattern.compile("^(\\d+);(\\d+);([^;]+);([0-9,]+)");

				BufferedReader reader = new BufferedReader(
						new InputStreamReader(new FileInputStream(previousBrokerFile)));
				try {
					reader.readLine(); // skip header row
					for (int lineCnt = 1; true; lineCnt++) {
						String line = reader.readLine();
						if (line == null) {
							break;
						}

						// Ignore incomplete histogram lines. A trailing comma indicates that the
						// histogram was not fully written, most likely because the previous run was
						// interrupted.
						Matcher m = lineRegexpPattern.matcher(line);
						if (m.matches() && !m.group(4).endsWith(",")) {
							int gainPlusLoss = Integer.parseInt(m.group(1));
							int iterations = Integer.parseInt(m.group(2));
							String filename = m.group(3);
							int[] histogram = parseStringFormOfIntarray(m.group(4));

							boolean ok1 = histogram.length == gainPlusLoss + 1;
							boolean ok2 = sum(histogram) == iterations;

							if (ok1 & ok2) {
								// registerResult(gainPlusLoss, histogram, iterations, filename);

								brokerFileWriter.write(gainPlusLoss + ";" + iterations + ";" + filename + ";"
										+ stringFormOf(histogram) + "\n");

								// Retrieve the task from the central task store. At this stage, all tasks are
								// still kept in the free-task map.
								BrokerTaskUnit u = taskMapFree.get(gainPlusLoss);
								assert (u != null);
								u.addHistogramToBrokerTaskUnit(histogram, iterations);
								// It is still in the map. I don't need to put
							} else {
								System.out.println("BAD LINE in Broker file (line: " + lineCnt + ") CODE=1");

							}
						} else {
							System.out.println("BAD LINE in Broker file (line: " + lineCnt + ") CODE=2");

						}

					} // for
				} finally {
					reader.close();
				}

				BrokerTaskUnit[] taskArray = taskMapFree.values().toArray(new BrokerTaskUnit[taskMapFree.size()]);
				taskMapFree.clear();
				for (int j = 0; j < taskArray.length; j++) {
					BrokerTaskUnit u = taskArray[j];
					u.updateInterestingGainIndex(maxIteration);

					int min = u.getMinimalInterestingCount();
					if (min > minIteration | u.cntTrial >= maxTrial) {
						// If all required gain counts have enough samples, or the trial limit has been
						// reached, mark the task as finished.
						taskListFinished.add(u);
					} else {
						taskMapFree.put(u.gainPlusLoss, u);
					}

				}

			} else {

				brokerFileWriter = new OutputStreamWriter(new FileOutputStream(brokerFile), "UTF-8");
				brokerFileWriter.write("gainPlusLoss;iterations;filename;histogramArray\n");
			}
		}

		@Override
		public String toString() {
			return "minIteration:" + minIteration + "  maxIteration:" + maxIteration + "\nFINISHED: " + taskListFinished
					+ "\n----------------\nWORKING:" + taskMapUnderWork + "\n------------\nFREE:" + taskMapFree;
		}

		private long timestampAtLastLog = 0;// System.currentTimeMillis();

		/*
		 * Write the log file only if at least five minutes have passed since the
		 * previous log update.
		 */
		private void tryToLog() {
			long t = System.currentTimeMillis();
			if (t - timestampAtLastLog > 5 * 60 * 1000) // Write a log update at most once every five minutes.
			{
				try {
					writeLogFile();
				} catch (IOException e) {
					e.printStackTrace();
				}

				timestampAtLastLog = t;
				// This method is only called from synchronized methods, so updating this field
				// is safe here.

			}

		}

		private void writeLogFile() throws FileNotFoundException, IOException {
			OutputStreamWriter writer = new OutputStreamWriter(new FileOutputStream(new File(logFileName)), "UTF-8");
			try {
				System.out.println("Writing task-broker log file: " + logFileName);
				String timeStamp = new SimpleDateFormat("yyyy.MM.dd.HH.mm.ss").format(new Date());
				writer.write(timeStamp + "\n");

				writer.write("minIteration:" + minIteration + "  maxIteration:" + maxIteration + "  maxTrial:"
						+ maxTrial + "\n");
				writer.write("FINISHED n=" + taskListFinished.size() + "+\n");
				for (BrokerTaskUnit btu : taskListFinished) {
					writer.write("  " + btu.toString(minIteration, maxIteration, maxTrial) + "\n");
				}
				writer.write("WORKING ON n=" + taskMapUnderWork.size() + "+\n");
				for (BrokerTaskUnit btu : taskMapUnderWork.values()) {
					writer.write("  " + btu.toString(minIteration, maxIteration, maxTrial) + "\n");
				}

				writer.write("FREE TASKS n=" + taskMapFree.size() + "+\n");
				for (BrokerTaskUnit btu : taskMapFree.values()) {
					writer.write("  " + btu.toString(minIteration, maxIteration, maxTrial) + "\n");
				}
			} finally {
				writer.close();
			}
		}


		// -------------------------------------------

		/**
		 * Register the result of a simulation task completed by a worker thread.
		 *
		 * A worker receives a task for a fixed gain+loss value, generates a block of
		 * random tree labelings, and writes the accepted labelings to a ZIP file. This
		 * method records the result in the broker file and updates the in-memory task
		 * state.
		 *
		 * @param gainPlusLoss total number of changes in the generated labelings
		 * @param histogram    histogram of generated labelings by gain count;
		 *                     histogram[g] stores how many generated labelings had g
		 *                     gains and gainPlusLoss - g losses
		 * @param iterations   number of random labelings attempted in this block
		 * @param filename     name of the ZIP file containing the accepted labelings,
		 *                     or null if no accepted labelings were written
		 */
		synchronized void registerResult(int gainPlusLoss, int[] histogram, int iterations, String filename)
				throws IOException {

			assert (gainPlusLoss + 1 == histogram.length); // The histogram is indexed by gain count. range: [0, 1, ...,
															// gainPlusLoss]
			assert (Arrays.stream(histogram).sum() == iterations); // TODO Verify whether mirrored gain/loss counts
																	// affect this assertion.

			// The ZIP file has already been written at this point. If the program is
			// interrupted externally, an unreferenced ZIP file may remain in the output
			// directory. However, the broker file should never reference a ZIP file that
			// has
			// not yet been written.

			brokerFileWriter
					.write(gainPlusLoss + ";" + iterations + ";" + filename + ";" + stringFormOf(histogram) + "\n");
			brokerFileWriter.flush();

			// Remove the task from the set of tasks currently under work.
			BrokerTaskUnit u = taskMapUnderWork.remove(gainPlusLoss);
			assert (u != null);

			u.addHistogramToBrokerTaskUnit(histogram, iterations);
			u.updateInterestingGainIndex(maxIteration);

			// Find the smallest number of accepted samples among the gain counts that are
			// still relevant for this task.
			int min = u.getMinimalInterestingCount();
			if (min > minIteration | u.cntTrial >= maxTrial) {
				// Mark the task as finished if enough samples have been collected for all
			    // relevant gain counts, or if the trial limit has been reached.
			 	taskListFinished.add(u);
			} else {
				// If the task remains active, keep only those gain counts for which we still
				// need additional samples.
				taskMapFree.put(gainPlusLoss, u);
			}

		}

		// -------------------------------------------

		synchronized MyTaskList.Row getTask() {
			tryToLog();

			// Take one task from the free-task map, move it to the under-work map, and
			// return it to the myWorker thread for processing.

			Iterator<Entry<Integer, BrokerTaskUnit>> iterator = taskMapFree.entrySet().iterator();
			if (iterator.hasNext()) {
				iterator.hasNext(); // TODO Check whether this no-op call is still needed; it may be a legacy iterator workaround.
				Entry<Integer, BrokerTaskUnit> mapEntry = iterator.next();
				BrokerTaskUnit u = taskMapFree.remove(mapEntry.getKey());
				taskMapUnderWork.put(mapEntry.getKey(), u);
				return new Row(u.gainPlusLoss, -1, u.interestingGainIdx);
			} else {
				return null;
			}
		}

		public synchronized void flush() throws IOException {
			System.out.println("Flushing the broker file.");
			brokerFileWriter.flush();
		}

		public synchronized void close() throws IOException {
			brokerFileWriter.close();
			System.out.println("Closing the broker file.");
		}

		/**
		 * @return comma-separated representation of the array values, e.g. "5,7,22,33"
		 */
		private static String stringFormOf(int[] arr) {
			StringBuilder sb = new StringBuilder(arr.length * 5);
			for (int i = 0; i < arr.length; i++) {
				sb.append(arr[i]).append(',');
			}

			if (arr.length > 0) {
				// Remove the final trailing comma.
				sb.setLength(sb.length() - 1);
			}
			return sb.toString();
		}

		/**
		 * Parse a comma-separated integer array representation, e.g. "5,7,22,33".
		 */
		private static int[] parseStringFormOfIntarray(String s) {
			String[] sa = s.split(",");
			int[] ia = new int[sa.length];
			for (int i = 0; i < sa.length; i++) {
				ia[i] = Integer.parseInt(sa[i]);
			}
			return ia;
		}

		static int sum(int[] arr) {
			int sum = 0; // initialize sum
			int i;

			for (i = 0; i < arr.length; i++)
				sum += arr[i];

			return sum;
		}

	} // end class TaskBroker

	/** return type */
	private static class ResultPack {
		int[] resultGainCntHistogram;
		HashMap<Integer, ArrayList<LabelingGandAndLoss>> resultMap;

		public ResultPack(int[] resultGainCntHistogram, HashMap<Integer, ArrayList<LabelingGandAndLoss>> resultMap) {
			super();
			this.resultGainCntHistogram = resultGainCntHistogram;
			this.resultMap = resultMap;
		}

	}

	/**
	 * Generate a block of random tree labelings for a fixed total number of
	 * changes.
	 *
	 * The returned histogram contains all generated labelings by gain count, while
	 * resultMap stores only those labelings whose gain count is marked as relevant
	 * in okGainCnt.
	 */
	static ResultPack findSomeGoodLabeling(int nGainPlusLoss, boolean[] okGainCnt, DoubleArrayTree tree, int iteration,
			Random random) {
		final int nn = nGainPlusLoss;

		int[] resultGainCntHistogram = new int[nn + 1];
		HashMap<Integer, ArrayList<LabelingGandAndLoss>> resultMap = new HashMap<Integer, ArrayList<LabelingGandAndLoss>>();

		IntegerStack gray = new IntegerStack(tree.size());
		// Reuse this stack in each iteration to avoid repeated memory allocation.

		for (int ii = 0; ii < iteration; ii++) {

			// Select a random set of nn tree nodes where state changes will occur.
			// Each selected node represents the child endpoint of an edge with either a
			// gain or a loss.

			// Set nn randomly selected positions to true.
			boolean[] idxOfChanges = generateRandomBooleanVector(tree.size(), nn, random);

			///////////////////////////////////////////
			// in this code block the tree will be traversed
			// and count how many gains and losses are induced by the
			// selected change positions (idxOfChaneges) , assuming that the gene is 1 at
			/////////////////////////////////////////// the root.
			int cntOfGains = 0;
			int cntOfLosses = 0;
			int[] geneValueMap = new int[tree.size()];

			gray.clean();
			gray.push(tree.getRoot());
			Arrays.fill(geneValueMap, -1);// TODO ez lehet, h felesleges
			geneValueMap[tree.getRoot()] = 1; // Assume that the gene is present at the root.

			while (!gray.isEmpty()) {
				int x = gray.pop();
				int parent = tree.getParent(x);
				int parentValue = parent == -1 ? 1 : geneValueMap[parent];
				boolean isChange = idxOfChanges[x];
				Integer value;
				if (isChange) {
					if (parentValue == 0) {
						value = 1;
						cntOfGains++;
					} else {
						value = 0;
						cntOfLosses++;
					}

				} else {
					value = parentValue;
				}

				// boolean notAlreadyContaind = (geneValueMap[x]== -1);
				geneValueMap[x] = value;

				gray.pushChildren(tree, x);
			} // while

			// Record the gain-count distribution of all generated labelings. The broker
			// uses this histogram to update progress, including labelings that were not
			// saved.
			resultGainCntHistogram[cntOfGains]++;

			if (okGainCnt[cntOfGains]) {
				// Traverse the tree again and record the node indices where gains and losses
				// occur.
				int[] gainArray = new int[cntOfGains];
				int[] lossArray = new int[cntOfLosses];
				int headGainArray = 0;
				int headLossArray = 0;

				gray.clean();
				gray.push(tree.getRoot());

				while (!gray.isEmpty()) {
					int x = gray.pop();
					boolean isChange = idxOfChanges[x];
					if (isChange) {
						if (geneValueMap[x] == 0) {
							// LOSS
							lossArray[headLossArray] = x;
							headLossArray++;

						} else {
							// GAIN
							gainArray[headGainArray] = x;
							headGainArray++;
						}

					}
					gray.pushChildren(tree, x);
				} // end of while - 2nd tree traverse block

				// Store the accepted labeling.
				LabelingGandAndLoss res = new LabelingGandAndLoss(gainArray, lossArray);
				ArrayList<LabelingGandAndLoss> list = resultMap.get(cntOfGains);
				if (list == null) {
					list = new ArrayList<LabelingGandAndLoss>();
					resultMap.put(cntOfGains, list);
				}
				list.add(res);

			} // if gainCnt is OK

		} // for(ii)

		return new ResultPack(resultGainCntHistogram, resultMap);
	}

	private static final boolean[] generateRandomBooleanVector(final int size, final int numberOfTrueValues,
			final Random random) {
		boolean[] idxOfChanges = new boolean[size]; // default values are false
		int cnt = numberOfTrueValues;
		while (cnt > 0) {

			int i = random.nextInt(idxOfChanges.length - 1) + 1; // Exclude the root node. That is why the +1 is needed

			if (!idxOfChanges[i]) {
				cnt--;
				idxOfChanges[i] = true;
			}
		}
		return idxOfChanges;
	}

	static class MyWorker implements Runnable {

		final TaskBroker taskBroker;
		final DoubleArrayTree tree;
		final int iterations;
		final Random random;
		final File outDir;
		Row taskRow;
		String id;

		public MyWorker(TaskBroker taskBroker, DoubleArrayTree tree, int iterations, Random random, File outDir,
				String workerId) {
			super();
			this.taskBroker = taskBroker;

			this.tree = tree;
			this.iterations = iterations;
			this.random = random;
			this.outDir = outDir;
			this.id = workerId;
		}

		@Override
		public void run() {

			// Each worker does the loop repeatedly on its own thread:
			// - requests a task from the broker,
			// - processes it,
			// - writes the accepted labelings to a ZIP file
			// - and registers the result with the broker.
			try {
				System.out.println("[" + id + "] Requesting a task from the broker.");
				while (true) {
					taskRow = taskBroker.getTask();
					System.out.println("[" + id + "] Task request completed.");
					if (taskRow == null) {
						System.out.println("[" + id + "] Worker thread terminating.");
						return;
						// If the broker has no more tasks, this worker terminates.
						// This is the only normal exit point from the worker loop.
					}

					System.out.println("[" + id + "]Starting [NN=" + taskRow.gainPlusLoss + "]");
					long t0 = System.nanoTime();

					// The main simulation step.
					ResultPack rp = findSomeGoodLabeling(taskRow.gainPlusLoss, taskRow.gains, tree, iterations, random);

					// The rp ResultPack belongs to one gain+loss value.
					// rp contains a The histogram which reports how many generated labelings had
					// each gain count.
					// sum(rt.resultGainCntHistogram)=itearions
					// The taskRow.gains is a boolean[] which specifies which gain counts should be
					// saved in rp.resultMap.

					long t1 = System.nanoTime();
					long timeInSec = (t1 - t0) / (1000000000);
					int timeInMin = (int) Math.floor(timeInSec / 60.0);
					long remainingSeconds = timeInSec - timeInMin * 60;
					System.out.println("[NN=" + taskRow.gainPlusLoss + "]- TIME:" + timeInSec + "sec = " + timeInMin
							+ "min+" + remainingSeconds + "sec");

					// rp.resultGainCntHistogram

					if (rp.resultMap.isEmpty()) {
						System.out.println("[" + id + "] No accepted labelings found. No ZIP file will be written.");
						System.out.println("[" + id + "] Registering the result with the broker.");
						taskBroker.registerResult(taskRow.gainPlusLoss, rp.resultGainCntHistogram, iterations, null);
						System.out.println("[" + id + "] Result registration completed.");

					} else {

						// Add a random suffix to make output filenames unique.
						String randomSuffix = randomChars(random, 8);
						try {
							String filename = "NN" + taskRow.gainPlusLoss + "-" + randomSuffix + ".ser.zip";
							saveToZipFile(rp.resultMap, outDir, filename);

							// If writing the ZIP file fails, execution skips this part and jumps to the
							// exception handler, where
							// an all-zero result is registered.
							System.out.println("[" + id + "] Registering the result with the broker.");

							taskBroker.registerResult(taskRow.gainPlusLoss, rp.resultGainCntHistogram, iterations,
									filename);
							System.out.println("[" + id + "] Result registration completed.");

						} catch (IOException e) {
							e.printStackTrace();
							System.out.println("[" + id + "] Registering an all-zero result with the broker.");
							taskBroker.registerResult(taskRow.gainPlusLoss, new int[rp.resultGainCntHistogram.length],
									0, null);
							System.out.println("[" + id + "] All-zero result registration completed.");
						}
					}
				} // while(true)
			} catch (IOException e) {
				e.printStackTrace();
			}

		}// end of run()

		@Override
		public String toString() {
			return "Worker_" + id + " working on " + taskRow;
		}

	} // end of class MyWorker

	/** Serialize resultMap into a ZIP file. */
	private static void saveToZipFile(Map<Integer, ArrayList<LabelingGandAndLoss>> resultMap, File dir, String filename)
			throws IOException, FileNotFoundException {

		File file = new File(dir, filename);
		OutputStream os = new FileOutputStream(file);
		ZipOutputStream zos = new ZipOutputStream(os);
		try {
			ZipEntry ze = new ZipEntry(filename);
			zos.putNextEntry(ze); // start first file in zip archive

			@SuppressWarnings("resource")
			ObjectOutputStream oos = new ObjectOutputStream(zos);

			oos.writeObject(resultMap);
			oos.flush();
			zos.closeEntry();
		} finally {
			zos.close();
		}
	}

	@SuppressWarnings("unchecked")
	public static Map<Integer, ArrayList<LabelingGandAndLoss>> loadFromZipFile(File file)
			throws IOException, FileNotFoundException {

		FileInputStream fileInputStream = new FileInputStream(file);
		ZipInputStream zipInputStream = new ZipInputStream(fileInputStream);
		try {
			// ZipEntry zipEntry =
			zipInputStream.getNextEntry();
			ObjectInputStream objectInputStream = new ObjectInputStream(zipInputStream);
			Map<Integer, ArrayList<LabelingGandAndLoss>> map;
			try {
				map = (Map<Integer, ArrayList<LabelingGandAndLoss>>) objectInputStream.readObject();
			} catch (ClassNotFoundException e) {
				throw new RuntimeException(e);
			}
			return map;

		} finally {
			zipInputStream.close();
		}
	}

	// ---------------------
	/**
	 * Generate a random character suffix for unique output filenames.
	 */
	private static String randomChars(Random random, int size) {
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < size; i++) {
			sb.append((char) ('a' + random.nextInt(24)));
		}
		return sb.toString();
	}

	////////////////////////////////////////////////////////
	public static void main(String[] args) throws IOException {
		System.out.println(Arrays.toString(args));
		if (args.length != 10) {
			System.err.println("Expected 10 command-line arguments:\n" + "  [1] numOfWorkers\n"
					+ "      Number of parallel worker threads.\n" + "  [2] treeNwkFilename\n"
					+ "      Path to the input phylogenetic tree in Newick format.\n" + "  [3] rootNodeName\n"
					+ "      Name of the root node of the subtree to analyse. Use 'root' for the full tree.\n"
					+ "  [4] iterationsPerBlock\n"
					+ "      Number of random tree labelings generated in one block for a given gain+loss value.\n"
					+ "  [5] minIterations\n"
					+ "      Minimum number of accepted samples to collect for each target gain/loss pair.\n"
					+ "  [6] maxIterations\n"
					+ "      Maximum number of samples to store for each target gain/loss pair. Additional hits\n"
					+ "      for that pair are ignored to save disk space.\n" + "  [7] maxTrial\n"
					+ "      Hard upper limit on the number of random labelings attempted for one task. The task\n"
					+ "      stops after this many trials even if minIterations has not been reached.\n"
					+ "  [8] taskFileName\n"
					+ "      Tab-separated task table produced by script-03-create_simulation_task_table.R.\n"
					+ "  [9] outDirName\n"
					+ "      Output directory for pregenerated simulation samples and the broker file.\n"
					+ "  [10] logFileName\n"
					+ "      Log file where the program periodically writes task-broker progress information.\n");
			System.exit(1);
		} else {
			int numOfWorkers = Integer.parseInt(args[0]);
			String treeNwkFilename = args[1].trim();
			String rootNodeName = args[2].trim();
			int iterationsPerBlock = Integer.parseInt(args[3]);
			int minIterations = Integer.parseInt(args[4]);
			int maxIterations = Integer.parseInt(args[5]);
			int maxTrial = Integer.parseInt(args[6]);//
			String taskFileName = args[7].trim();
			String outDirName = args[8].trim();
			String logFileName = args[9].trim();

			File brokerFile = new File(outDirName, "_BROKER_FILE.txt");

			long randomInitializationCounter = 398347539;

			File taskFile = new File(taskFileName);
			Util.createParentDirectory(taskFile);

			ArrayList<Row> taskList = MyTaskList.load(taskFile);
			// System.out.println(taskList);
			System.out.println("taskList loaded.");

			System.out.println("OutDir: " + outDirName);
			File outDir = new File(outDirName);
			outDir.mkdirs();

			TaskBroker taskBroker = new TaskBroker(taskList, minIterations, maxIterations, maxTrial, logFileName,
					brokerFile);

			taskBroker.tryToLog();

			DoubleArrayTree tree;
			if ("root".equalsIgnoreCase(rootNodeName)) {
				tree = new DoubleArrayTree(new File(treeNwkFilename));

			} else {
				tree = new DoubleArrayTree(new File(treeNwkFilename), rootNodeName);
			}

			long randomSeedBase = System.nanoTime();

			Thread[] workerThreadArray = new Thread[numOfWorkers];
			for (int i = 0; i < numOfWorkers; i++) {

				randomInitializationCounter--;
				Random random = new Random(randomSeedBase + randomInitializationCounter);

				workerThreadArray[i] = new Thread(
						new MyWorker(taskBroker, tree, iterationsPerBlock, random, outDir, "Worker_" + i),
						"Worker_" + i);
			}

			for (int i = 0; i < numOfWorkers; i++) {
				workerThreadArray[i].start();
			}

			for (int i = 0; i < numOfWorkers; i++) {
				workerThreadArray[i].setPriority(Thread.MIN_PRIORITY);
			}

			// Wait for all worker threads to finish.
			for (int i = 0; i < numOfWorkers; i++) {

				try {
					workerThreadArray[i].join();
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
			}

			taskBroker.writeLogFile();

			taskBroker.close();
			System.out.println("Done.");
		}
	}// main() method

}
