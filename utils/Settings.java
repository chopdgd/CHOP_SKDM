package utils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.PrintStream;
import java.util.Iterator;
import java.util.TreeMap;
import java.util.Map.Entry;

public class Settings {
	public static final String VERSION = "beta";
	
	public static String WORKING_DIR = "";
	
	protected static int DEC_PLACES_PERC = 2;
	protected static int DEC_PLACES_PVAL = 5;
	protected static double PVAL_CUTOFF = .05;
	public static boolean RESIDUE_DISCOVERY_BY_DELTA = false;
	public static int POSITIONS_TO_INTERROGATE = 0;
	public static boolean PRINT_SIGNIFICANT_ONLY = false;
	public static String MANUAL_STRATIFICATION = "";
	
	public static final double GET_PVAL_CUTOFF(){
		return PVAL_CUTOFF;
	}
	public static final int GET_DEC_PLACES_PERC(){
		return DEC_PLACES_PERC;
	}
	public static final int GET_DEC_PLACES_PVAL(){
		return DEC_PLACES_PVAL;
	}
	public static final int GET_POSITIONS_TO_INTERROGATE(){
		return POSITIONS_TO_INTERROGATE;
	}
	public static final boolean GET_PRINT_SIGNIFICANT_ONLY(){
		return PRINT_SIGNIFICANT_ONLY;
	}
	public static final String GET_MANUAL_STRATIFICATION(){
		return MANUAL_STRATIFICATION;
	}
	
	public static final boolean PRINT_ON = true;
	public static final String ARROW = " ***";
	public static final String BLANK = "-------";

    //changed : to = since alleles now have the ":" in their ids.
	public static final String EMPTY_STRING = " = = ";
	
	public static final boolean SUPPRESS_WARNINGS = false;
	
	private static final String[] DEFAULT_LABELS = {"a: DECIMAL PLACES FOR PERCENTAGES",
													"b: DECIMAL PLACES FOR P-VALUES",
													"c: P-VALUE THRESHOLD",
													"d: BOOST P-VALUE METHOD",
													"e: POSITIONS TO INTERROGATE",
													"f: PRINT SIGNIFICANT ONLY",
													"g: MANUAL STRATIFICATION"
	};
	private static final String[] DEFAULT_TOOLTIPS = {"number of decimals to output for percentages",
													  "number of decimals to output for p-values",
													  "significance threshold",
													  "run an alternative analysis (based on deltas) that results in a lower correction",
													  "interrogate this many positions from the start of the alignment [0 = all]",
													  "hide items whose corrected p-value is less than the threshold",
													  "to be implemented..."
	};
	private static final String[] DEFAULT_VALUES = {"2",
			   										"5",
			   										".05",
			   										"false",
			   										"0",
			   										"true",
			   										""
	};
	
	public static final void load(){
		File fpath = new File(".");
		String path;
		try {
			WORKING_DIR = fpath.getCanonicalPath() + File.separatorChar;
		} catch (Exception e) {}
		path = WORKING_DIR + "settings.txt";
		
		try{
			BufferedReader br;
			String inputLine, thisSetting;
			File f = new File(path);
			int line = 0;
			if (f.canRead()){
				br = new BufferedReader(new FileReader(new File(path)));
	        	while ((inputLine = br.readLine()) != null && inputLine.length() > 1){
	        		thisSetting = inputLine.split(":")[2].trim();
	        		switch(line){
	        		case 0:
	        			DEC_PLACES_PERC = Integer.parseInt(thisSetting);
	        			if (DEC_PLACES_PERC < 0 || DEC_PLACES_PERC > 10)
	        				DEC_PLACES_PERC = Integer.parseInt(DEFAULT_VALUES[line]);
	        			break;
	        		case 1:
	        			DEC_PLACES_PVAL = Integer.parseInt(thisSetting);
	        			if (DEC_PLACES_PVAL < 0 || DEC_PLACES_PVAL > 10)
	        				DEC_PLACES_PVAL = Integer.parseInt(DEFAULT_VALUES[line]);
	        			break;
	        		case 2:
	        			PVAL_CUTOFF = Double.parseDouble(thisSetting);
	        			if (PVAL_CUTOFF < 0 || PVAL_CUTOFF > 1)
	        				PVAL_CUTOFF = Integer.parseInt(DEFAULT_VALUES[line]);
	        			break;
	        		case 3:
	        			RESIDUE_DISCOVERY_BY_DELTA = thisSetting.equalsIgnoreCase("true")?true:false;
	        			break;
	        		case 4:
	        			POSITIONS_TO_INTERROGATE = Integer.parseInt(thisSetting);
	        			if (POSITIONS_TO_INTERROGATE<0)
	        				POSITIONS_TO_INTERROGATE = Integer.parseInt(DEFAULT_VALUES[line]);
	        			break;
	        		case 5:
	        			PRINT_SIGNIFICANT_ONLY = thisSetting.equalsIgnoreCase("true")?true:false;
	        			break;
	        		case 6:
	        			MANUAL_STRATIFICATION = thisSetting;
	        			break;
	        		}
	        		line++;
	        	}
	        	br.close();
			}
        }
        catch (Exception e){
        	Tools.print("\nError in settings file. Reverting to defaults...\n" + e.getMessage(), "font", "color=red");
        	setSettings(null);
        }
	}
	
	public static final void setSettings(Iterator it){
		String strSettings = "";
		try{
			if (it == null){
				for (int i = 0; i < DEFAULT_VALUES.length; i++){
					strSettings += DEFAULT_LABELS[i] + " (" + DEFAULT_TOOLTIPS[i] + ") : " + DEFAULT_VALUES[i] + "\n";
				}
			}
			else{
				int i = 0;
				String thisSetting;
				while (it.hasNext()){
					thisSetting = it.next().toString();
					strSettings += DEFAULT_LABELS[i] + " (" + DEFAULT_TOOLTIPS[i] + ") : " + thisSetting + "\n";
					i++;
				}
			}
			FileOutputStream fout;		
		    fout = new FileOutputStream(WORKING_DIR + "settings.txt");
		    new PrintStream(fout).print(strSettings);
		    fout.close();
		}catch (Exception e){
			setSettings(null);
		}
	}
	
	public static final Iterator<Entry<String,Object>> getSettings(){
		TreeMap<String, Object> ret = new TreeMap<String, Object>();
		int i = 0;
		ret.put(DEFAULT_LABELS[i] + "!" + DEFAULT_TOOLTIPS[i++], DEC_PLACES_PERC);
		ret.put(DEFAULT_LABELS[i] + "!" + DEFAULT_TOOLTIPS[i++], DEC_PLACES_PVAL);
		ret.put(DEFAULT_LABELS[i] + "!" + DEFAULT_TOOLTIPS[i++], PVAL_CUTOFF);
		ret.put(DEFAULT_LABELS[i] + "!" + DEFAULT_TOOLTIPS[i++], RESIDUE_DISCOVERY_BY_DELTA);
		ret.put(DEFAULT_LABELS[i] + "!" + DEFAULT_TOOLTIPS[i++], POSITIONS_TO_INTERROGATE);
		ret.put(DEFAULT_LABELS[i] + "!" + DEFAULT_TOOLTIPS[i++], PRINT_SIGNIFICANT_ONLY);
		ret.put(DEFAULT_LABELS[i] + "!" + DEFAULT_TOOLTIPS[i++], MANUAL_STRATIFICATION);
		return ret.entrySet().iterator();
	}
}
