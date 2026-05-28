package ger.maddison;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;


public class MyTaskList
{
    public static class Row
    {
	public final int gainPlusLoss;
	public final int n;
	public final boolean[] gains;

	public Row(int gain_plus_loss, int n, boolean[] gains)
	{
	    super();
	    this.gainPlusLoss = gain_plus_loss;
	    this.n = n;
	    this.gains = gains;
	}

	@Override
	public String toString()
	{
	    StringBuilder sb = new StringBuilder();
	    sb.append("[NN=").append(gainPlusLoss);
	    sb.append("; n=").append(n);
	    sb.append("; gains=");
	    for (int i = 0; i < gains.length; i++)
	    {
		if (gains[i])
		{
		    sb.append(i).append(',');
		}
	    }
	    sb.setLength(sb.length() - 1);
	    sb.append("]\n");
	    return sb.toString();
	}
    }

    public static ArrayList<Row> load(File file) throws IOException
    {

	ArrayList<Row> list = new ArrayList<Row>();
	BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(file)));
	try
	{
	    in.readLine();// skip header
	    String line;
	    while (true)
	    {
		line = in.readLine();
		if (line == null)
		{
		    return list;
		}
		String[] f = line.split("\t");
		int gain_plus_loss = Integer.parseInt(f[0]);
		int n = Integer.parseInt(f[1]);
		boolean[] gains = new boolean[gain_plus_loss + 1];
		String[] ff = f[2].split(",");
		for (String xs : ff)
		{
		    int x = Integer.parseInt(xs);
		    gains[x] = true;
		}
		list.add(new Row(gain_plus_loss, n, gains));

	    }

	} finally
	{
	    in.close();
	}
    }
}
