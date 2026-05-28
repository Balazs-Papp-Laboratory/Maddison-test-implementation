package ger.maddison;

import java.io.File;

public class Util
{
    public static void createParentDirectory(File f)
    {
	File dir = f.getParentFile();

	if (dir != null)
	{
	    if (!dir.exists())
	    {
		dir.mkdirs();
	    }
	}
    }
}
