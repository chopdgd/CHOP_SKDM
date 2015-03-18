package main;

import utils.FisherExact;
import utils.Settings;
import utils.Tools;
import web.HLAAlign;

import java.util.*;
import java.util.Map.Entry;

public class Residue {
	private HLAAlign hlaalign;
	private ArrayList<String> ctrl;
	private ArrayList<String> pat;
	private int ctrlCount;
	private int patCount;
	private HashMap<String, ResidueEntry> pvalMap;
	private TreeMap<String, ResidueEntry> sortedMap;
	private int comparisons, pocketComps;
	private Integer[] pocketList;
	
	private static final boolean DETAILS = true;
	
	private boolean PRINT_ALLELES = !Settings.GET_PRINT_SIGNIFICANT_ONLY();
	
	@SuppressWarnings("unchecked")
	public Residue(HLAAlign hlaalign, ArrayList<String> ctrl, int ctrlCount, ArrayList<String> pat, int patCount){
		this.hlaalign = hlaalign;
		this.ctrl = ctrl;
		this.pat = pat;
		this.ctrlCount = ctrlCount;
		this.patCount = patCount;
		pvalMap = new HashMap<String, ResidueEntry>();
		loadPocketArray();
		checkRes();
		sortedMap = new TreeMap<String, ResidueEntry>(new ResidueEntryComparator());
		sortedMap.putAll(pvalMap);
	}
	
	private final void checkRes(){
		ArrayList<String> alleles = hlaalign.getAlleles();
		ArrayList<Double> pvals = hlaalign.getPvals();
		ArrayList<Double> ors = hlaalign.getORs();
		if (alleles.size() < 1) return;
		ArrayList<String> tempAlleles = new ArrayList<String>(alleles.size() - 1);
		char[][] aligns = hlaalign.getAligns();
		char[] refarr = hlaalign.getRefSeq();
		int i = 0, j = 0, break_point = aligns.length;
		comparisons = 0; pocketComps = 0;
		int pok = 0;
		int interrogateUpTo = (Settings.GET_POSITIONS_TO_INTERROGATE() == 0? aligns[0].length : Math.min(Settings.GET_POSITIONS_TO_INTERROGATE(), aligns[0].length));
		if (DETAILS)
			Tools.printDetails("\n");
		
		// get negative breakpoint
		if (Settings.RESIDUE_DISCOVERY_BY_DELTA) break_point = hlaalign.getNegativeBreak();
		if (break_point < 0) return;
		
		// get reference sequence
		int refi = hlaalign.getRefSeqIndex();
		if (refi >= 0){
			// replace reference sequence with ----...
			for (i = 0; i < aligns[refi].length; i++)
				aligns[refi][i] = '-';
		}
		
		// start scan
		char curr;
		String aasf, aasr;
		
		// FORWARD
		for (j = 0; j < interrogateUpTo; j++){
			aasf = ""; aasr = "";
			for (i = 0; i < break_point; i++){
				curr = aligns[i][j];
				// if AA not collected, collect
				if (aasf.indexOf(curr) == -1 && curr != '*'){
					aasf += curr;
				}
			}
			// if in negatives, uncollect
			for (i = break_point; i < aligns.length; i++){
				curr = aligns[i][j];
				if (aasf.indexOf(curr) > -1){
					aasf = aasf.replaceAll(String.valueOf(curr), "");
				}
			}
			// if we're left with AAs... ============================================
			if (Settings.RESIDUE_DISCOVERY_BY_DELTA? aasf.length() > 0 : aasf.length() > 1){
				boolean association;
				for (int aacount = 0; aacount < aasf.length(); aacount++){
					association = true; // positively associated
					
					// collect alleles
					for (i = 0; i < break_point; i++){
						curr = aligns[i][j];
						if (curr == aasf.charAt(aacount)){
							tempAlleles.add(alleles.get(i)); // collect alleles
						}
					}
					// replace dash with AA name from reference seq
					if (aasf.charAt(aacount) == '-') // do not fuckin touch
						aasf = aasf.replaceAll("-", String.valueOf(refarr[j]));
					String alsID = "";
					double pval = -1, or = -1;
					
					// if only 1, we know pval
					if (tempAlleles.size() == 1){
						int pvalIndex = alleles.indexOf(tempAlleles.get(0));
						pval = pvals.get(pvalIndex);
						or = ors.get(pvalIndex);
						if (hlaalign.getNegativeBreak() > 0 && pvalIndex>=hlaalign.getNegativeBreak())
							association = false; // negatively associated
						alsID = tempAlleles.get(0) + "|";
						if (!pvalMap.containsKey(alsID))
							pvalMap.put(alsID, new ResidueEntry(new String[] {tempAlleles.get(0)}, pval, or, association, j, aasf.charAt(aacount)));
						else
							pvalMap.put(alsID, pvalMap.get(alsID).add(j, aasf.charAt(aacount)));
						if (DETAILS){
							Tools.printDetails("\n#" + (comparisons+1) + "\t");
							Tools.printDetails(hlaalign.getLocus() + "_" + aasf.charAt(aacount) + "" + (j+1) + (association?"+":"-"), "b");
							Tools.printDetails("\t{see " + pvalMap.get(alsID).getResAlleles()[0] + " allele}");

                            //changed : to = in the statement below.
							Tools.printDetails("\t" + "P= " + Tools.formatDec(pval, Settings.GET_DEC_PLACES_PVAL()+1) + "\tOR= " + Tools.formatDec(or, Settings.GET_DEC_PLACES_PERC()+1));
						}
					}
					
					else{
						for (i = 0; i < tempAlleles.size(); i++){
							alsID += tempAlleles.get(i) + "|"; // build Allele ID
						}
						
						// if ID not contained, compute pval
						if (!pvalMap.containsKey(alsID)){
							int cn = 0, pn = 0;
							// controls
							for (int k = 0; k < ctrl.size(); k=k+2){
								i = 0;
								while (i < tempAlleles.size() && !tempAlleles.get(i).equals(ctrl.get(k)))
									i++;
								if (i < tempAlleles.size() && tempAlleles.get(i).equals(ctrl.get(k))){
									cn++;
								}
								else{
									i = 0;
									while (i < tempAlleles.size() && !tempAlleles.get(i).equals(ctrl.get(k+1)))
										i++;
									if (i < tempAlleles.size() && tempAlleles.get(i).equals(ctrl.get(k+1))){
										cn++;
									}
								}
								
							}
							// patients
							for (int k = 0; k < pat.size(); k=k+2){
								i = 0;
								while (i < tempAlleles.size() && !tempAlleles.get(i).equals(pat.get(k)))
									i++;
								if (i < tempAlleles.size() && tempAlleles.get(i).equals(pat.get(k))){
									pn++;
								}
								else{
									i = 0;
									while (i < tempAlleles.size() && !tempAlleles.get(i).equals(pat.get(k+1)))
										i++;
									if (i < tempAlleles.size() && tempAlleles.get(i).equals(pat.get(k+1))){
										pn++;
									}
								}
							}
							pval = FisherExact.exact22total(cn, ctrlCount, pn, patCount);
							or = Tools.oddsRatio(pn, patCount-pn, cn, ctrlCount-cn);
							if (pn*10000/patCount - cn*10000/ctrlCount < 0) association = false; // negatively associated
							
							if (DETAILS){
								Tools.printDetails("\n#" + (comparisons+1) + "\t");
								Tools.printDetails(hlaalign.getLocus() + "_" + aasf.charAt(aacount) + "" + (j+1) + (association?"+":"-"), "b");
								Tools.printDetails("\t{" + pn + ", " + (patCount-pn) + ", " + cn + ", " + (ctrlCount-cn) + "}");

                                //changed : to = in the statement below.
								Tools.printDetails("\t" + "P= " + Tools.formatDec(pval, Settings.GET_DEC_PLACES_PVAL()+1) + "\tOR= " + Tools.formatDec(or, Settings.GET_DEC_PLACES_PERC()+1));
							}
							
							String[] tempAllArr = new String[tempAlleles.size()];
							System.arraycopy(tempAlleles.toArray(), 0, tempAllArr, 0, tempAlleles.size());
							pvalMap.put(alsID, new ResidueEntry(tempAllArr, pval, or, association, j, aasf.charAt(aacount)));
						}
						else{
							pvalMap.put(alsID, pvalMap.get(alsID).add(j, aasf.charAt(aacount)));
							if (DETAILS){
								Tools.printDetails("\n#" + (comparisons+1) + "\t");
								Tools.printDetails(hlaalign.getLocus() + "_" + aasf.charAt(aacount) + "" + (j+1) + (association?"+":"-"), "b");
								Tools.printDetails("\t{same as " + pvalMap.get(alsID).getResAAs()[0] + "" + (pvalMap.get(alsID).getResPoss()[0]+1) + "}");

                                //changed : to = in the statement below.
								Tools.printDetails("\t" + "P= " + Tools.formatDec(pvalMap.get(alsID).getPval(), Settings.GET_DEC_PLACES_PVAL()+1) + "\tOR= " + Tools.formatDec(pvalMap.get(alsID).getOR(), Settings.GET_DEC_PLACES_PERC()+1));
							}
						}
					}
					comparisons++;
					tempAlleles.clear();
					if (pok < pocketList.length && (j+1) == pocketList[pok])
						pocketComps++;
				}
			}
			
			if (Settings.RESIDUE_DISCOVERY_BY_DELTA){
				// REVERSE
				for (i = aligns.length-1; i >= break_point; i--){
					curr = aligns[i][j];
					if (aasr.indexOf(curr) == -1 && curr != '*'){
						aasr += curr;
					}
				}
				for (i = break_point-1; i >= 0; i--){
					curr = aligns[i][j];
					if (aasr.indexOf(curr) > -1){
						aasr = aasr.replaceAll(String.valueOf(curr), "");
					}
				}
				// if we're left with AAs... ============================================
				if (aasr.length() > 0){
					for (int aacount = 0; aacount < aasr.length(); aacount++){
						// collect alleles
						for (i = aligns.length-1; i >= break_point; i--){
							curr = aligns[i][j];
							if (aasr.indexOf(curr) > -1){
								tempAlleles.add(alleles.get(i)); // collect alleles
							}
						}
						// replace dash with AA name from reference seq
						if (aasr.charAt(aacount) == '-')
							aasr = aasr.replaceAll("-", String.valueOf(refarr[j]));
						String alsID = "";
						double pval = -1, or = -1;
						
						// if only 1, we know pval
						if (tempAlleles.size() == 1){
							int pvalIndex = alleles.indexOf(tempAlleles.get(0));
							pval = pvals.get(pvalIndex);
							or = ors.get(pvalIndex);
							alsID = tempAlleles.get(0);
							if (!pvalMap.containsKey(alsID))
								pvalMap.put(alsID, new ResidueEntry(new String[] {tempAlleles.get(0)}, pval, or, false, j, aasr.charAt(aacount)));
							else
								pvalMap.put(alsID, pvalMap.get(alsID).add(j, aasr.charAt(aacount)));
							if (DETAILS){
								Tools.printDetails("\n#" + (comparisons+1) + "\t");
								Tools.printDetails(hlaalign.getLocus() + "_" + aasr.charAt(aacount) + "" + (j+1) + "-", "b");
								Tools.printDetails("\t{see " + pvalMap.get(alsID).getResAlleles()[0] + " allele}");

                                //changed : to = in the statement below.
								Tools.printDetails("\t" + "P= " + Tools.formatDec(pval, Settings.GET_DEC_PLACES_PVAL()+1) + "\tOR= " + Tools.formatDec(or, Settings.GET_DEC_PLACES_PERC()+1));
							}
						}
						
						// else compute!
						else{
							for (i = 0; i < tempAlleles.size(); i++){
								alsID += tempAlleles.get(i) + "|"; // build Allele ID
							}
							// if ID not contained, compute pval
							if (!pvalMap.containsKey(alsID)){
								int cn = 0, pn = 0;
								// controls
								for (int k = 0; k < ctrl.size(); k=k+2){
									i = 0;
									while (i < tempAlleles.size() &&
											!tempAlleles.get(i).equals(ctrl.get(k)))
										i++;
									if (i < tempAlleles.size() && tempAlleles.get(i).equals(ctrl.get(k)))
										cn++;
									else{
										i = 0;
										while (i < tempAlleles.size() && 
												!tempAlleles.get(i).equals(ctrl.get(k+1)))
											i++;
										if (i < tempAlleles.size() && tempAlleles.get(i).equals(ctrl.get(k+1)))
											cn++;
									}
								}
								// patients
								for (int k = 0; k < pat.size(); k=k+2){
									i = 0;
									while (i < tempAlleles.size() &&
											!tempAlleles.get(i).equals(pat.get(k)))
										i++;
									if (i < tempAlleles.size() && tempAlleles.get(i).equals(pat.get(k)))
										pn++;
									else{
										i = 0;
										while (i < tempAlleles.size() && 
												!tempAlleles.get(i).equals(pat.get(k+1)))
											i++;
										if (i < tempAlleles.size() && tempAlleles.get(i).equals(pat.get(k+1)))
											pn++;
									}
								}
								pval = FisherExact.exact22total(cn, ctrlCount, pn, patCount);
								or = Tools.oddsRatio(pn, patCount-pn, cn, ctrlCount-cn);
								
								if (DETAILS){
									Tools.printDetails("\n#" + (comparisons+1) + "\t");
									Tools.printDetails(hlaalign.getLocus() + "_" + aasr.charAt(aacount) + "" + (j+1) + "-", "b");
									Tools.printDetails("\t{" + pn + ", " + (patCount-pn) + ", " + cn + ", " + (ctrlCount-cn) +"}");

                                    //changed : to = in the statement below.
									Tools.printDetails("\t" + "P= " + Tools.formatDec(pval, Settings.GET_DEC_PLACES_PVAL()+1) + "\tOR= " + Tools.formatDec(or, Settings.GET_DEC_PLACES_PERC()+1));
								}
								
								String[] tempAllArr = new String[tempAlleles.size()];
								System.arraycopy(tempAlleles.toArray(), 0, tempAllArr, 0, tempAlleles.size());
								pvalMap.put(alsID, new ResidueEntry(tempAllArr, pval, or, false, j, aasr.charAt(aacount)));
							}
							else{
								pvalMap.put(alsID, pvalMap.get(alsID).add(j, aasr.charAt(aacount)));
								if (DETAILS){
									Tools.printDetails("\n#" + (comparisons+1) + "\t");
									Tools.printDetails(hlaalign.getLocus() + "_" + aasr.charAt(aacount) + "" + (j+1) + "-", "b");
									Tools.printDetails("\t{same as " + pvalMap.get(alsID).getResAAs()[0] + "" + (pvalMap.get(alsID).getResPoss()[0]+1) + "}");

                                    //changed : to = in the statement below.
									Tools.printDetails("\t" + "P= " + Tools.formatDec(pvalMap.get(alsID).getPval(), Settings.GET_DEC_PLACES_PVAL()+1) + "\tOR= " + Tools.formatDec(pvalMap.get(alsID).getOR(), Settings.GET_DEC_PLACES_PERC()+1));
								}
							}
						}
						comparisons++;
						tempAlleles.clear();
						if (pok < pocketList.length && (j+1) == pocketList[pok])
							pocketComps++;
					}
				}
			} // if RESIDUE_DISCOVERY_BY_DELTA
			if (pok < pocketList.length && (j+1) >= pocketList[pok])
				pok++;
		} // for j (columns)
	}
	
	private final void loadPocketArray(){
		//counting the DISCREET pocket positions!
		HashMap<Integer, Integer> hm = new HashMap<Integer, Integer>();
		int[][] pockets = hlaalign.getLocus().getPockets();
		for (int n = 0; n < pockets.length; n++)
			for (int m = 0; m < pockets[n].length; m++)
				hm.put(pockets[n][m], 1);
		//sort
		Iterator <Entry<Integer, Integer>>it = hm.entrySet().iterator();
		TreeMap<Integer, Integer> tm = new TreeMap<Integer, Integer>();
		while (it.hasNext()){
			tm.put(it.next().getKey(), 1);
		}
		//to array
		int n = hm.size();
		Integer[] pp = new Integer[n];
		pp = tm.keySet().toArray(pp);
		pocketList = pp;
	}
	
	private final String checkPockets(){
		int[][] pockets = hlaalign.getLocus().getPockets();
		if (pockets.length < 1) return Settings.EMPTY_STRING;
		
		if (pocketList.length < 1) return Settings.EMPTY_STRING;
		int[] sigRes = new int[pockets.length];
		for(int i = 0; i < sigRes.length; i++) sigRes[i] = 0;
		StringBuffer[] sba = new StringBuffer[pockets.length];
		StringBuffer ret = new StringBuffer();

        //changed : to = in the statement below.

        //also changed "Tools.underline" from "=" to "_"
		ret.append("\n\n\n" + Tools.underline("=HLA-" + hlaalign.getLocus() + " Pocket Residues="));
		ret.append("\n> p-value correction is " + pocketComps + " (equal to pocket AAs interrogated)");
		for (int b = 0; b < pockets.length; b++){
			sba[b] = new StringBuffer();
			if (pockets[b].length > 0){
                System.out.println("\nPocket " + (b+1) + " [");
				sba[b].append("\nPocket " + (b+1) + " [");
				for (int p = 0; p < pockets[b].length; p++){
                    System.out.println(pockets[b][p]+(p<pockets[b].length-1?",":""));
					sba[b].append(pockets[b][p]+(p<pockets[b].length-1?",":""));
                }
				sba[b].append("]");
			}
		}
		int currPos = -1;
		Iterator<ResidueEntry> it = sortedMap.values().iterator();
		ResidueEntry curr;
		while (it.hasNext()){
			curr = it.next();
			// for each pos in residue entry block
			for (int j = 0; j < curr.getResPoss().length; j++){
				currPos = curr.getResPoss()[j] + 1;
				// run through pocket matrix:
				for (int p = 0; p < pockets.length; p++){
					for (int q = 0; q < pockets[p].length; q++){
						if (currPos == pockets[p][q]){
							if (Settings.PRINT_SIGNIFICANT_ONLY? curr.getPval()*pocketComps < Settings.GET_PVAL_CUTOFF() : curr.getPval() < Settings.GET_PVAL_CUTOFF()){
								if (sigRes[p]==0){
                                    System.out.println("\nPos\tAA\tAssoc\tp-val\tp^corr\tOR\n-------\t-------\t-------\t-------\t-------\t-------");
									sba[p].append("\nPos\tAA\tAssoc\tp-val\tp^corr\tOR\n-------\t-------\t-------\t-------\t-------\t-------");
                                }
                                System.out.println(curr.toString(pocketComps, false, j));
								sba[p].append(curr.toString(pocketComps, false, j));
								sigRes[p]++;
							}
						}
					}
				}
			}
			//lastPos = curr.getResPoss()[0] + 1;
		}
		
		boolean isEmpty = true;
		for (int b = 0; b < sba.length; b++){
			if (sba[b].length()>0){
				ret.append(sba[b].toString());
			}
			if (sigRes[b]>0){
				isEmpty = false;
			}
		}
		if (isEmpty){
			ret.append("\nNo significant residues found in HLA-" + hlaalign.getLocus() + " pockets!");
        }

		return ret.toString();
	}
	
	public final ResidueEntry[] getResidues(){
		return sortedMap.values().toArray(new ResidueEntry[sortedMap.size()]);
	}
	
	public final LinkedList<ResidueEntry> getResidueColl(){
		return new LinkedList<ResidueEntry>(sortedMap.values());
	}
	
	public final Locus getLocus(){
		return hlaalign.getLocus();
	}
	                    
	@SuppressWarnings("unchecked")
	public final String toString(){
		String loc = hlaalign.getLocus().toString();
		StringBuffer sb = new StringBuffer();
		double currPval = 1;

        //changed : to = in the statement below.
		sb.append("\n\nHLA-" + loc + " Residues=");
		sb.append("\n> p-value correction is " + comparisons + " (equal to number of AA interrogated)");
		Iterator<ResidueEntry> it = sortedMap.values().iterator();
		int i = 0;
		ResidueEntry curr;
		boolean isEmpty = true;
		while (it.hasNext()){
			curr = it.next();
			currPval = curr.getPval();
			// apply correction?
			//currPval = Settings.RESIDUE_DISCOVERY_PVAL_CORRECTION? currPval * comparisons : currPval;


			if (Settings.PRINT_SIGNIFICANT_ONLY? currPval * (double)comparisons < Settings.GET_PVAL_CUTOFF() : currPval < Settings.GET_PVAL_CUTOFF()){
				if (isEmpty)
					sb.append("\n" + (PRINT_ALLELES?"Alls\t":"")
                            + "Pos\t\tAA\tAssoc\tp-val\tp^corr\tOR\n" + (PRINT_ALLELES?"-------\t":"") + "--------------\t-------\t-------\t-------\t-------\t-------");
				sb.append(curr.toString(comparisons, PRINT_ALLELES));
				isEmpty = false;
			}
			i++;
		}
		if (isEmpty)
			sb.append("\nNo significant residues found for HLA-" + hlaalign.getLocus() + "!");
		return sb.toString() + checkPockets();
	}
	
	public final class ResidueEntryComparator implements Comparator {
		public final int compare(Object resID1, Object resID2) {
			int ret;
			ResidueEntry r1 = pvalMap.get((String)resID1);
			ResidueEntry r2 = pvalMap.get((String)resID2);
			int a1 = r1.getResPoss()[0];
			int a2 = r2.getResPoss()[0];
			if (a1>=a2)
				ret = 1;
			else ret = -1;
			return ret;
		}
	}
}