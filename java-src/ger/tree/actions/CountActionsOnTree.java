package ger.tree.actions;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.Map;

import ger.maddison.IntMatrix;
import ger.maddison.Node;
import ger.maddison.NodeUtils;
import ger.maddison.NwkNode;

/**
 * Count binary gene presence/absence states and gain/loss transitions on a
 * tree.
 *
 * This command-line program reads a phylogenetic tree in Newick format and a
 * DAT-ordered binary gene presence/absence matrix produced by the previous
 * preprocessing step of the pipeline. The input matrix contains one row per
 * gene or orthogroup and one column per tree node. Matrix values are expected
 * to be binary:
 *
 * 0 = gene absent 
 * 1 = gene present
 *
 * The program iterates over all genes in the matrix and counts, for each gene,
 * how many tree edges belong to each of the four possible parent-to-child
 * state-transition categories:
 *
 * 0 -> 0 absent 
 * 0 -> 1 gain 
 * 1 -> 0 loss 
 * 1 -> 1 present
 *
 * The analysis can be restricted to a subtree by providing the name of the root
 * node of that subtree. Use "root" to analyse the full tree.
 *
 * Command-line arguments:
 *
 * [1] path to the Newick tree file [2] name of the root node of the subtree to
 * analyse; use "root" for the full tree [3] path to the DAT-ordered binary gene
 * presence/absence matrix [4] output TSV filename
 *
 * The output is a tab-separated table with one row per gene and the following
 * columns:
 *
 * gene1 absent gain loss present
 *
 * The column name "gene1" is kept for compatibility with downstream pipeline
 * steps, although it represents the gene or orthogroup identifier.
 */
public class CountActionsOnTree {
	public static void main(String[] args) throws IOException {
		if (args.length != 4) {
			System.err.println("Expected 4 command-line arguments:\n" + //
					"  [1] path to the Newick tree file\n" + //
					"  [2] name of the root node of the subtree to analyse; use 'root' for the full tree\n" + //
					"  [3] path to the DAT-ordered binary gene presence/absence matrix\n" + //
					"  [4] output TSV filename\n");

			if (args.length != 0) {
				System.err.println("Received " + args.length + " command-line arguments:");
				for (int i = 0; i < args.length; i++) {
					System.err.println("  \"" + args[i] + "\"\n");
				}
			}

			System.exit(1);

		} else {
			String treeFilePath = args[0].trim();
			String rootNodeName = args[1].trim();
			String presenceMatrixFilePath = args[2].trim();
			String outFileName = args[3].trim();
			// -------------

			System.out.println("Loading tree.");
			Node tree = NwkNode.parseNWK(new File(treeFilePath));
			if (!"root".equalsIgnoreCase(rootNodeName)) {
				tree = NodeUtils.findFirst(tree, rootNodeName);
			}
			System.out.println("done.");
			int treeSize = NodeUtils.collectAllDescendants(tree).size();
			System.out.println("tree size=" + treeSize + ";  root name=" + tree.getName());

			System.out.println("Loading matrix.");
			IntMatrix x = new IntMatrix();
			x.initFromTsv(new File(presenceMatrixFilePath));

			System.out.println("Data loaded");
			// ---------------

			File outFile = new File(outFileName);
			File dir = outFile.getParentFile();
			if (dir != null) {
				dir.mkdirs();
			}
			writeTransitionCounts(outFile, tree, x);

			System.out.println("Done");

		}
	}

	private static void writeTransitionCounts(File outFile, Node tree, IntMatrix x)
			throws IOException, FileNotFoundException {

		long time0 = System.currentTimeMillis();

		BufferedWriter out = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outFile), "UTF-8"));
		try {
			out.write("gene1\tabsent\tgain\tloss\tpresent\n"); // header row of the out file

			for (int i0 = 0; i0 < x.nrow; i0++) {

				if (i0 % 4000 == 100) {
					long time1 = System.currentTimeMillis();
					double elapsed = time1 - time0;
					double timeLeft = elapsed * ((double) (x.nrow - i0) / i0); // milisecounds

					System.out.println(outFile.getName() + " TimeLeft: " + Math.round(timeLeft / 60000) + " min ");

				}

				Map<String, Integer> map1 = x.getRowAsMap(i0);

				out.write(x.rownames[i0]);
				// out.write('\t');

				int[] res = NodeUtils.countLossAndGainEvents(map1, tree);
				for (int j1 = 0; j1 < res.length; j1++) {
					out.write('\t');
					out.write(String.valueOf(res[j1]));
				}
				out.write('\n');

			}
		} finally {
			out.close();
		}
	}
}
