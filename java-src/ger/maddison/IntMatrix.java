package ger.maddison;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Map;

public class IntMatrix
{
    public String[] colnames;
    public String[] rownames;
    public int[][] values;
    public int ncol;
    public int nrow;

    public void initFromTsv(File tsvFile) throws IOException
    {
	initFromSeparatedTextFile(tsvFile, '\t');
    }

    public void initFromCsv(File csvFile) throws IOException
    {
	initFromSeparatedTextFile(csvFile, ',');
    }

    public void initFromSeparatedTextFile(File file, char separator) throws IOException
    {

	BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(file)));
	try
	{
	    String headerLine = in.readLine();
	    int t1 = headerLine.indexOf(separator);
	    colnames = headerLine.substring(t1 + 1).split("" + separator);
	    ncol = colnames.length;

	    int cnt = 0;
	    while (true)
	    {
		String line = in.readLine();
		if (line == null)
		{
		    break;
		}
		cnt++;
	    }
	    nrow = cnt;
	    rownames = new String[nrow];
	    values = new int[nrow][ncol];

	} finally

	{
	    in.close();
	}

	int lineNumber = -1;
	in = new BufferedReader(new InputStreamReader(new FileInputStream(file)));
	try
	{

	    in.readLine();
	    for (lineNumber = 0; lineNumber < nrow; lineNumber++)
	    {
		String line = in.readLine();
		int t2 = line.indexOf(separator);
		rownames[lineNumber] = line.substring(0, t2);
		String[] tmp = line.substring(t2 + 1).split("" + separator);
		for (int j = 0; j < ncol; j++)
		{
		    values[lineNumber][j] = Integer.parseInt(tmp[j]);
		}
	    }

	} catch (RuntimeException e)
	{
	    System.err.println("Error in line " + (lineNumber + 1));
	    throw e;
	}

	finally

	{
	    in.close();
	}

    }

    public Map<String, Integer> getRowAsMap(int idx)
    {
	HashMap<String, Integer> map = new HashMap<String, Integer>(ncol);
	int[] row = values[idx];
	for (int i = 0; i < ncol; i++)
	{
	    map.put(colnames[i], row[i]);
	}

	return map;
    }

    public int[] getRow(int idx)
    {
	int[] row = values[idx];

	return row;
    }

    
    public static void main(String[] args) throws IOException
    {
	IntMatrix x = new IntMatrix();
	x.initFromTsv(new File("data/RES_ALL_count_gain2_result_FAMILY.tsv"));
	System.out.println(x);
    }
}
