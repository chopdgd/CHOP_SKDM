package main;

import utils.Settings;
import utils.Tools;

import java.util.*;

public class ListUnique {	
	private LinkedList<String> alleles = new LinkedList<String>();
	private HashMap<String, Double[]> unique;
	private Locus loc;
	private int population;
	private ArrayList<String> sortedAlleles;
	private String title;
	
	@SuppressWarnings("unchecked")
	public ListUnique(Locus loc, LinkedList<String> alleles, int population, String title){
		this.alleles = alleles;
		this.loc = loc;
		this.population = population;
		this.title = title;
		calcUnique();
		TreeMap<String, Double[]> sorted = new TreeMap<String, Double[]>(new AlleleMapComparator(unique, 1, 1));
		sorted.putAll(unique);
		sortedAlleles = new ArrayList<String>(sorted.keySet());
	}
	
	private void calcUnique(){
		unique = new HashMap<String, Double[]>();
		String s = null, ns = null;
		boolean flag = false, end = false;
		Iterator<String> it = alleles.iterator();
		while (!end){
			int popCount = 0, popPerc = 0, allCount = 0, allPerc = 0;
			if (!flag)
				if (it.hasNext())
					s = it.next();
				else{
					s = "";
					end = true;
				}
			else{
				s = ns;
				flag = false;
				if (!it.hasNext())
					end = true;
			}
			if (s.length() > 0){
				if (unique.containsKey(s)){
					popCount = unique.get(s)[0].intValue();
					allCount = unique.get(s)[2].intValue();
				}
				if (it.hasNext()){
					if ((ns = it.next()).equals(""))
						allCount++;
					else
						flag = true;
				}
				popCount++;
				allCount++;
				popPerc = (int)Math.round(popCount*100*Math.pow(10,Settings.GET_DEC_PLACES_PERC())/population);
				allPerc = (int)Math.round(allCount*100*Math.pow(10,Settings.GET_DEC_PLACES_PERC())/(population*2));
				Double[] data = {(double)popCount, (double)popPerc, (double)allCount, (double)allPerc};
				unique.put(s, data);
			}
		}
	}
		
	public String[] getAlleles(){
		return (String[])unique.keySet().toArray(new String[unique.size()]);
	}
	
	public HashMap<String, Double[]> getAllelesMap(){
		return unique;
	}
	
	public ArrayList<String> getAllelesList(){
		return new ArrayList<String>(unique.keySet());
	}
	
	public String[] getOriginalAlleles(){
		return alleles.toArray(new String[alleles.size()]);
	}
	
	/*public String[] getAllelesSorted(){
		return sortedAlleles.toArray(new String[sortedAlleles.size()]);
	}*/
	
	public int getPopulation(){
		return population;
	}
	
	public String getTitle(){
		return title;
	}
	
	public String toString(){
		StringBuffer sb = new StringBuffer();

        ////changed : to = in the statement below.
		sb.append("\n\n" + (loc == Locus.ERR?"":"HLA-") + loc + " summary=\nAllele\tPop Freq\tAllele Freq");
		String popPerc, allPerc;
		for (int i = 0; i < sortedAlleles.size(); i++){
			popPerc = Tools.formatDec(unique.get(sortedAlleles.get(i))[1].intValue(), Settings.GET_DEC_PLACES_PERC());
			allPerc = Tools.formatDec(unique.get(sortedAlleles.get(i))[3].intValue(), Settings.GET_DEC_PLACES_PERC());
			sb.append("\n" + sortedAlleles.get(i) + "\t" + unique.get(sortedAlleles.get(i))[0].intValue() + "  " + popPerc + "%\t" + unique.get(sortedAlleles.get(i))[2].intValue() + "  " + allPerc + "%");
		}
		sb.append("\n" + sortedAlleles.size() + "\tunique alleles total.\n" + population + "\tsamples total.");
		return sb.toString();
	}
}