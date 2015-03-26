package src.main;

import java.util.Iterator;
import java.util.LinkedList;

import src.utils.Settings;
import src.utils.Tools;

public class ResidueEntry {

	private String[] resAlleles;
	private double pval;
	private double or;
	private boolean positiveAssociation;
	private LinkedList<Integer> residuePos;
	private LinkedList<Character> residueAA;
	//TODO: private int[] twobytwo;  //store twobytwo and compute everything here
	
	public ResidueEntry(String[] resAlleles, double pval, double or, boolean isPositive, int resPos, char resAA){
		this.resAlleles = resAlleles;
		this.pval = pval;
		this.or = or;
		positiveAssociation = isPositive;
		residuePos = new LinkedList<Integer>();
		residuePos.add(resPos);
		residueAA = new LinkedList<Character>();
		residueAA.add(resAA);
	}

	public final String[] getResAlleles(){
		return resAlleles;
	}
	
	public final boolean isPositive(){
		return positiveAssociation;
	}

	public final double getPval(){
		return pval;
	}
	
	public final double getOR(){
		return or;
	}

	public final int[] getResPoss(){
		Iterator<Integer> it = residuePos.iterator();
		int[] ret = new int[residuePos.size()];
		int i = 0;
		while (it.hasNext())
			ret[i++] = it.next();
		return ret;
	}
	
	public final String printResPoss(){
		StringBuffer sb = new StringBuffer();
		int[] rp = getResPoss();
		for (int i = 0; i < rp.length; i++)
			sb.append((rp[i]+1)+(i==rp.length-1?"":","));
		return sb.toString();
	}
	
	public final char[] getResAAs(){
		Iterator<Character> it = residueAA.iterator();
		char[] ret = new char[residueAA.size()];
		int i = 0;
		while (it.hasNext())
			ret[i++] = it.next();
		return ret;
	}
	public final String printResAAs(){
		StringBuffer sb = new StringBuffer();
		char[] rp = getResAAs();
		for (int i = 0; i < rp.length; i++)
			sb.append(rp[i]);
		return sb.toString();
	}
	
	public final ResidueEntry add(int pos, char res){
		residuePos.add(pos);
		residueAA.add(res);
		return this;
	}
	
	public final String toString(){
		return toString(1, true);
	}
	
	public final String toString(int pValCorrection, boolean printAlleles){
		StringBuffer sb = new StringBuffer();
		double currPval = getPval();
		double currPvalCorr = (currPval * (double)pValCorrection > 1? (double)1 : currPval * (double)pValCorrection);
		double currOR = getOR();
		if (currPval > Math.pow(10, Settings.GET_DEC_PLACES_PVAL())) currPval = (int)Math.pow(10, Settings.GET_DEC_PLACES_PVAL());
		if (printAlleles){
			for (int n = 0; n < getResAlleles().length; n++){
				sb.append("\n" + getResAlleles()[n]);
				if (n == 0)
					sb.append("\t" + (printResPoss().length()<8?printResPoss()+"\t":printResPoss()) + "\t" + printResAAs() + "\t" + (isPositive()? "+" : "-") + "\t" + Tools.formatDec(currPval, Settings.GET_DEC_PLACES_PVAL()) + "\t" + Tools.formatDec(currPvalCorr, Settings.GET_DEC_PLACES_PVAL()) + "\t" + Tools.formatDec(currOR, Settings.GET_DEC_PLACES_PERC()) + (currPvalCorr < Settings.GET_PVAL_CUTOFF()? Settings.ARROW:""));
			}
			sb.append("\n" + Settings.BLANK + "\t" + Settings.BLANK + Settings.BLANK + "\t" + Settings.BLANK + "\t" + Settings.BLANK + "\t" + Settings.BLANK + "\t" + Settings.BLANK + "\t" + Settings.BLANK);
		}
		else{
			sb.append("\n" + (printResPoss().length()<8?printResPoss()+"\t":printResPoss()) + "\t" + printResAAs() + "\t" + (isPositive()? "+" : "-") + "\t" + Tools.formatDec(currPval, Settings.GET_DEC_PLACES_PVAL()) + "\t" + Tools.formatDec(currPvalCorr, Settings.GET_DEC_PLACES_PVAL()) + "\t" + Tools.formatDec(currOR, Settings.GET_DEC_PLACES_PERC()) + (currPvalCorr < Settings.GET_PVAL_CUTOFF()? Settings.ARROW:""));
		}
		return sb.toString();
	}
	
	public final String toString(int pValCorrection, boolean printAlleles, int index){
		StringBuffer sb = new StringBuffer();
		double currPval = getPval();
		double currPvalCorr = (currPval * (double)pValCorrection > 1? (double)1 : currPval * (double)pValCorrection);
		double currOR = getOR();
		if (currPval > Math.pow(10, Settings.GET_DEC_PLACES_PVAL())) currPval = (int)Math.pow(10, Settings.GET_DEC_PLACES_PVAL());
		if (printAlleles){
			for (int n = 0; n < getResAlleles().length; n++){
				sb.append("\n" + getResAlleles()[n]);
				if (n == 0)
					sb.append("\t" + (getResPoss()[index]+1) + "\t" + getResAAs()[index] + "\t" + (isPositive()? "+" : "-") + "\t" + Tools.formatDec(currPval, Settings.GET_DEC_PLACES_PVAL()) + "\t" + Tools.formatDec(currPvalCorr, Settings.GET_DEC_PLACES_PVAL()) + "\t" + Tools.formatDec(currOR, Settings.GET_DEC_PLACES_PERC()) + (currPvalCorr < Settings.GET_PVAL_CUTOFF()? Settings.ARROW:""));
			}
			sb.append("\n" + Settings.BLANK + "\t" + Settings.BLANK + "\t" + Settings.BLANK + "\t" + Settings.BLANK + "\t" + Settings.BLANK + "\t" + Settings.BLANK + "\t" + Settings.BLANK);
		}
		else{
			sb.append("\n" + (getResPoss()[index]+1) + "\t" + getResAAs()[index] + "\t" + (isPositive()? "+" : "-") + "\t" + Tools.formatDec(currPval, Settings.GET_DEC_PLACES_PVAL()) + "\t" + Tools.formatDec(currPvalCorr, Settings.GET_DEC_PLACES_PVAL()) + "\t" + Tools.formatDec(currOR, Settings.GET_DEC_PLACES_PERC()) + (currPvalCorr < Settings.GET_PVAL_CUTOFF()? Settings.ARROW:""));
		}
		return sb.toString();
	}
}
