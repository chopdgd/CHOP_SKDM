package main;

import utils.FisherExact;
import utils.Settings;
import utils.Tools;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.TreeMap;

public class Delta {
	private ListUnique ctrl, pat;
	private HashMap <String, Double[]>del;
	public Locus loc;
	private ArrayList<String> sortedAlleles;
	
	@SuppressWarnings("unchecked")
	public Delta(ListUnique ctrl, ListUnique pat, Locus loc){
		this.ctrl = ctrl;
		this.pat = pat;
		this.loc = loc;
		del = new HashMap<String, Double[]>();
		deltas();
		TreeMap<String, Double[]> sorted = new TreeMap<String, Double[]>(new AlleleMapComparator(del, 0, 1));
		sorted.putAll(del);
		sortedAlleles = new ArrayList<String>(sorted.keySet());
	}
	
	private void deltas(){
		HashMap <String, Double[]>c = ctrl.getAllelesMap();
		HashMap <String, Double[]>p = pat.getAllelesMap();
		String[] cval = ctrl.getAlleles();
		ArrayList <String>pval = pat.getAllelesList();
		Double[] data; //percentage delta, p-value
		for (int i = 0; i < cval.length; i++){
			if (p.containsKey(cval[i])){
				data = new Double[]{(double)(p.get(cval[i])[1] - c.get(cval[i])[1]), FisherExact.exact22total(c.get(cval[i])[0], ctrl.getPopulation(), p.get(cval[i])[0], pat.getPopulation()), Tools.oddsRatio(p.get(cval[i])[0], pat.getPopulation()-p.get(cval[i])[0], c.get(cval[i])[0], ctrl.getPopulation()-c.get(cval[i])[0])};
				del.put(cval[i], data);
				pval.remove(cval[i]);
			}
			else{
				data = new Double[]{(double)(- c.get(cval[i])[1]), FisherExact.exact22total(c.get(cval[i])[0], ctrl.getPopulation(), 0, pat.getPopulation()), Tools.oddsRatio(0, pat.getPopulation(), c.get(cval[i])[0], ctrl.getPopulation()-c.get(cval[i])[0])};
				del.put(cval[i], data);
			}
		}
		for (int i = 0; i < pval.size(); i++){
			data = new Double[]{(double)p.get(pval.get(i))[1], FisherExact.exact22total(0, ctrl.getPopulation(), p.get(pval.get(i))[0], pat.getPopulation()), Tools.oddsRatio(p.get(pval.get(i))[0], pat.getPopulation()-p.get(pval.get(i))[0], 0, ctrl.getPopulation())};
			del.put(pval.get(i), data);
		}
		//System.out.println("ctrl: " + ctrl.getPopulation() + "\t pat: " + pat.getPopulation());
	}
	
	public ArrayList<String> getAllelesSorted(){
		return sortedAlleles;
	}
	
	public ArrayList<String> getDeltas(){
		ArrayList<String> ret = new ArrayList<String>(sortedAlleles.size());
		for (int i = 0; i < sortedAlleles.size(); i++){
			ret.add(i, Tools.formatDec(del.get(sortedAlleles.get(i))[0].intValue(), Settings.GET_DEC_PLACES_PERC()));
		}
		return ret;
	}
	
	public ArrayList<Double> getPvals(){
		ArrayList<Double> ret = new ArrayList<Double>(sortedAlleles.size());
		for (int i = 0; i < sortedAlleles.size(); i++){
			ret.add(i,del.get(sortedAlleles.get(i))[1]);
		}
		return ret;
	}
	
	public ArrayList<Double> getORs(){
		ArrayList<Double> ret = new ArrayList<Double>(sortedAlleles.size());
		for (int i = 0; i < sortedAlleles.size(); i++){
			ret.add(i,del.get(sortedAlleles.get(i))[2]);
		}
		return ret;
	}
	
	public String toString(){
		StringBuffer sb = new StringBuffer();
        //changed : to = in the statement below.
		sb.append("\n\nDelta between " + pat.getTitle() + " and " + ctrl.getTitle() + " for locus HLA-" + loc + "=\nAllele\tDelta\tp^corr\tOR");
		Iterator<String> it = sortedAlleles.iterator();
		int allelesTotal = sortedAlleles.size();
		while (it.hasNext()){
			String curr = it.next();
			sb.append("\n" + curr + "\t" + Tools.formatDec(del.get(curr)[0].intValue(), Settings.GET_DEC_PLACES_PERC()) + "%\t" + Tools.formatDec(Math.min(del.get(curr)[1]*(double)allelesTotal,1), Settings.GET_DEC_PLACES_PVAL()) + "\t" + Tools.formatDec(del.get(curr)[2], Settings.GET_DEC_PLACES_PERC()) + (del.get(curr)[1]*(double)allelesTotal < Settings.GET_PVAL_CUTOFF()? Settings.ARROW : ""));
		}
		sb.append("\n" + allelesTotal + "\talleles total.");
		return sb.toString();
	}
}