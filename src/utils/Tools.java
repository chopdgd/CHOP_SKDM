package src.utils;

import src.gui.Window;
import src.main.Run;

import java.text.DecimalFormat;

public class Tools {
	
	public static final String formatDec(int n, int DEC_PLACES){
		String perc = String.valueOf(n);
		boolean neg = false;
		if (perc.length()>1 && perc.charAt(0) == '-'){
			neg = true;
			perc = perc.substring(1);
		}
		if (DEC_PLACES < 1) return perc;
		for (int i = perc.length(); i <= DEC_PLACES; i++)
			perc = "0" + perc;
		return (neg?"-":"") + perc.substring(0, perc.length()-DEC_PLACES) + "." + perc.substring(perc.length()-DEC_PLACES);
	}
	
	public static final String formatDec(double n, int DEC_PLACES){
		DEC_PLACES = Math.min(DEC_PLACES, String.valueOf(Double.MAX_VALUE).length()-2);
		if (n == (double)0)
			return "0";
		String template = "0.";
		for (int i = 0; i < DEC_PLACES; i++)
			template += "#";
		if (n < (double)0.001)
			template = template.substring(0, (template.length()-2>3? template.length()-2 : template.length())) + "E0";
		DecimalFormat df = new DecimalFormat(template);
		return df.format(n);
	}
	
	public static final double oddsRatio(double n11, double n12, double n21, double n22){
		double ahalf = (double)0.5;
		double fwd = (n11+ahalf) * (n22+ahalf) / ((n12+ahalf) * (n21+ahalf));
		//double rev = (n12+ahalf) * (n21+ahalf) / ((n11+ahalf) * (n22+ahalf));
		return fwd;//Math.max(fwd, rev);
	}
	
	public static final double oddsRatio(int n11, int n12, int n21, int n22){
		return oddsRatio((double)n11, (double)n12, (double)n21, (double)n22);
	}
	
	public static final double oddsRatio(double[] darr){
		return oddsRatio(darr[0], darr[1], darr[2], darr[3]);
	}
	
	public static final double oddsRatio(int[] narr){
		return oddsRatio((double)narr[0], (double)narr[1], (double)narr[2], (double)narr[3]);
	}
	
	public static final String underline(String s, char uline){
		String ret = s + "\n";
		for (int i = 0; i < s.length(); i++)
			ret += uline;
		return ret;
	}
	
	public static final String underline(String s){
		return underline(s, '_');
	}

	public static final void printAll(String s){
		if (Settings.PRINT_ON){
			System.out.print(s);
			Window.appendOutputTxt(s.replaceAll("\n", "<br>"));
			Window.appendDetailsTxt(s.replaceAll("\n", "<br>"));
		}
		Run.LOG.append(s);
	}
	
	public static final void printAll(String s, String tag, String attribs){
		if (Settings.PRINT_ON){
			System.out.print(s);
			Window.appendOutputTxt("<" + tag + (attribs.length()>0?" " + attribs:"") + ">" + s.replaceAll("\n", "<br>") + "</" + tag + ">");
			Window.appendDetailsTxt("<" + tag + (attribs.length()>0?" " + attribs:"") + ">" + s.replaceAll("\n", "<br>") + "</" + tag + ">");
		}
		Run.LOG.append(s);
	}
	
	public static final void print(String s){
		if (Settings.PRINT_ON){
			System.out.print(s);
			Window.appendOutputTxt(s.replaceAll("\n", "<br>"));
		}
		Run.LOG.append(s);
	}
	
	public static final void print(String s, String tag, String attribs){
		if (Settings.PRINT_ON){
			System.out.print(s);
			Window.appendOutputTxt("<" + tag + (attribs.length()>0?" " + attribs:"") + ">" + s.replaceAll("\n", "<br>") + "</" + tag + ">");
		}
		Run.LOG.append(s);
	}
	
	public static final void printErr(String s){
		if (Settings.PRINT_ON){
			System.out.print(s);
			Window.appendOutputTxt("<font color=red><b>" + s.replaceAll("\n", "<br>") + "</b></font>");
		}
		Run.LOG.append(s);
	}
	
	public static final void print(String s, String tag){
		print(s, tag, "");
	}
	
	public static final void printDetails(String s){
		if (Settings.PRINT_ON){
			Window.appendDetailsTxt(s.replaceAll("\n", "<br>"));
		}
		Run.LOG.append(s);
	}
	
	public static final void printDetails(String s, String tag, String attribs){
		if (Settings.PRINT_ON){
			Window.appendDetailsTxt("<" + tag + (attribs.length()>0?" " + attribs:"") + ">" + s.replaceAll("\n", "<br>") + "</" + tag + ">");
		}
		Run.LOG.append(s);
	}
	
	public static final void printDetails(String s, String tag){
		printDetails(s, tag, "");
	}
	
	public static final String pop(String s){
		return s.substring(0, s.length()-1);
	}
	
	public static final String pop(String s, int n){
		if (n>0 && n < s.length())
			return s.substring(0, s.length()-n);
		else
			return pop(s);
	}
}
