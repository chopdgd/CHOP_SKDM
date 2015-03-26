package src.main;

import src.utils.FisherExact;
import src.utils.Settings;
import src.utils.Tools;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.TreeMap;

public class LinkageDisequilibrium {

	private ListUnique[] pats;
	private ListUnique[] ctrls;
	private Residue[] residues;
	private int[][] popMatrix; //built along the way, lists subjects VS present/absent AA (1,0) or 2 if present twice (homozygote)
	private ArrayList<int[][]> alleleMatrix;
	private int[] patToCtrl;
	private int popPatToCtrl;
	private int[] hlaToNext;
	private LinkedList<LDEntry> patCtrlEntries;
	private LinkedList<String> zygosities;
	private ResidueEntry[][] resArr;
	private int whichMatrix;
	
	private static final double tests = 5;
	private static final double additionalTests = 8;
	private static final double zygosityTests = 3;
	private static final int MISSING = -9;
	
	private static final boolean DETAILS = true;
	private static final boolean DEBUG = false;
	
	private static final boolean ALLOW_POS_NEG_INTERACTION = true;
	
	public LinkageDisequilibrium(ListUnique[] pats, ListUnique[] ctrls, Residue[] residues){
		this.pats = pats;
		this.ctrls = ctrls;
		this.residues = residues;
		patCtrlEntries = new LinkedList<LDEntry>();
		popMatrix = new int[0][0];
		patToCtrl = new int[residues.length];
		hlaToNext = new int[residues.length+1];
		hlaToNext[0] = 0;
		alleleMatrix = new ArrayList<int[][]>(residues.length);
		zygosities = new LinkedList<String>();
		buildMatrices();
	}
	
	private void buildMatrices(){
		resArr = new ResidueEntry[residues.length][];
		for (int hla = 0; hla < residues.length; hla++){
			int[][] patCtrlMatrix;
			int[][] patCtrlPopMatrix;
			
			// remove insignificant residues
			LinkedList<ResidueEntry> resList = residues[hla].getResidueColl();
			LinkedList<ResidueEntry> resListClean = new LinkedList<ResidueEntry>();
			Iterator<ResidueEntry> it = resList.iterator();
			while (it.hasNext()){
				ResidueEntry re = it.next();
				if (re.getPval() < Settings.GET_PVAL_CUTOFF())
					resListClean.add(re);
			}
			
			resArr[hla] = new ResidueEntry[resListClean.size()];
			resArr[hla] = resListClean.toArray(resArr[hla]);

			
			// init
			String[] patAlleles = pats[hla].getAlleles();
			String[] ctrlAlleles = ctrls[hla].getAlleles();
			String[] patPopAlleles = pats[hla].getOriginalAlleles();
			String[] ctrlPopAlleles = ctrls[hla].getOriginalAlleles();
			
			patCtrlMatrix = new int[resArr[hla].length][patAlleles.length+ctrlAlleles.length];
			patCtrlPopMatrix = new int[resArr[hla].length][patPopAlleles.length/2+ctrlPopAlleles.length/2];
			
			patToCtrl[hla] = patAlleles.length;
			popPatToCtrl = patPopAlleles.length/2;
			
			for (int i = 0; i < patCtrlMatrix.length; i++)
				for (int j = 0; j < patCtrlMatrix[0].length; j++)
					patCtrlMatrix[i][j] = 0;
			for (int i = 0; i < patCtrlPopMatrix.length; i++)
				for (int j = 0; j < patCtrlPopMatrix[0].length; j++)
					patCtrlPopMatrix[i][j] = 0;
			//////////////////////
			
			// fill
			// PATIENTS //////////////////////
			for (int i = 0; i < resArr[hla].length; i++){
				String[] currAlleles = resArr[hla][i].getResAlleles();
				for (int j = 0; j < patAlleles.length; j++){
					int k = 0;
					while (k < currAlleles.length && !currAlleles[k].equals(patAlleles[j])){
						k++;
					}
					if (k < currAlleles.length){ // we have a hit
						patCtrlMatrix[i][j] = 1;
					}
				}
			}
			for (int i = 0; i < resArr[hla].length; i++){
				String[] currAlleles = resArr[hla][i].getResAlleles();
				for (int j = 0; j < patPopAlleles.length; j=j+2){
					if (patPopAlleles[j].equals("")){
						patCtrlPopMatrix[i][j/2] = MISSING; // missing data
					}
					else{
						int z = 0;
						for (int k = 0; k < currAlleles.length; k++){
							if (currAlleles[k].equals(patPopAlleles[j]))
								z++;
							if (currAlleles[k].equals(patPopAlleles[j+1]))
								z++;
						}
						patCtrlPopMatrix[i][j/2] = z;
					}
				}
			}
			
			// CONTROLS //////////////////////
			for (int i = 0; i < resArr[hla].length; i++){
				String[] currAlleles = resArr[hla][i].getResAlleles();
				for (int j = 0; j < ctrlAlleles.length; j++){
					int k = 0;
					while (k < currAlleles.length && !currAlleles[k].equals(ctrlAlleles[j])){
						k++;
					}
					if (k < currAlleles.length){ // we have a hit
						patCtrlMatrix[i][j+patToCtrl[hla]] = 1;
					}
				}
			}
			for (int i = 0; i < resArr[hla].length; i++){
				String[] currAlleles = resArr[hla][i].getResAlleles();
				for (int j = 0; j < ctrlPopAlleles.length; j=j+2){
					if (ctrlPopAlleles[j].equals("")){
						patCtrlPopMatrix[i][j/2+popPatToCtrl] = MISSING; // missing data
					}
					else{
						int z = 0;
						for (int k = 0; k < currAlleles.length; k++){
							if (currAlleles[k].equals(ctrlPopAlleles[j]))
								z++;
							if (currAlleles[k].equals(ctrlPopAlleles[j+1]))
								z++;
						}
						patCtrlPopMatrix[i][j/2+popPatToCtrl] = z;
					}
				}
			}
			hlaToNext[hla+1] = hlaToNext[hla] + patCtrlPopMatrix.length;
			//////////////////////
			
			// collect matrices
			alleleMatrix.add(patCtrlMatrix);
			addToPopMatrix(patCtrlPopMatrix);
		}
		
		// TODO: Manual Stratification!
		if (!Settings.GET_MANUAL_STRATIFICATION().equals("")){
			
		}
		
		// print the matrices?
		if (DEBUG){
			int e = 0;
			for (int m = 0; m < alleleMatrix.size(); m++){
				Tools.printDetails("\n\nPatient Matrix (" + (m+1) + ")\n");
				for (int i = 0; i < alleleMatrix.get(m).length; i++){
					for (int j = 0; j < patToCtrl[m]; j++){
						if (i == 0 && j == 0) {
							Tools.printDetails("\t");
							for (int k = 0; k < patToCtrl[m]; k++){
								Tools.printDetails(pats[m].getAlleles()[k] + " ");
							}
							Tools.printDetails("\ncounts\t");
							for (int k = 0; k < patToCtrl[m]; k++){
								Tools.printDetails(String.format("%1$-5s", pats[m].getAllelesMap().get(pats[m].getAlleles()[k])[2].intValue()));
							}
							Tools.printDetails("\n");
						}
						if (j == 0)
							Tools.printDetails(resArr[m][i].getResAAs()[0] + "" + (resArr[m][i].getResPoss()[0]+1) + "\t");
						Tools.printDetails(alleleMatrix.get(m)[i][j] + "    ");
					}
					Tools.printDetails("\n");
				}
				Tools.printDetails("\nControl Matrix (" + (m+1) + ")\n");
				for (int i = 0; i < alleleMatrix.get(m).length; i++){
					for (int j = patToCtrl[m]; j < alleleMatrix.get(m)[0].length; j++){
						if (i == 0 && j == patToCtrl[m]) {
							Tools.printDetails("\t");
							for (int k = patToCtrl[m]; k < alleleMatrix.get(m)[0].length; k++){
								Tools.printDetails(ctrls[m].getAlleles()[k-patToCtrl[m]] + " ");
							}
							Tools.printDetails("\ncounts\t");
							for (int k = patToCtrl[m]; k < alleleMatrix.get(m)[0].length; k++){
								Tools.printDetails(String.format("%1$-5s", ctrls[m].getAllelesMap().get(ctrls[m].getAlleles()[k-patToCtrl[m]])[2].intValue()));
							}
							Tools.printDetails("\n");
						}
						if (j == patToCtrl[m])
							Tools.printDetails(resArr[m][i].getResAAs()[0] + "" + (resArr[m][i].getResPoss()[0]+1) + "\t");
						Tools.printDetails(alleleMatrix.get(m)[i][j] + "    ");
					}
					Tools.printDetails("\n");
				}
			}
			e = 0;
			Tools.printDetails("\n\nPopulation Joint Matrix [pat | ctrl]\n");
			for (int i = 0; i < popMatrix.length; i++){
				for (int j = 0; j < popMatrix[0].length; j++){
					if (i == 0 && j == 0){
						Tools.printDetails("\t");
						for (int k = 0; k < popPatToCtrl; k++){
							Tools.printDetails(String.format("%1$-4s", (k+1)));
						}
						Tools.printDetails(" |  ");
						for (int k = popPatToCtrl; k < popMatrix[0].length; k++){
							Tools.printDetails(String.format("%1$-4s", (k-popPatToCtrl+1)));
						}
						Tools.printDetails("\n");
					}
					if (j == 0){
						while (i >= hlaToNext[e+1]){
							e++;
						}
						Tools.printDetails(resArr[e][i-hlaToNext[e]].getResAAs()[0] + "" + (resArr[e][i-hlaToNext[e]].getResPoss()[0]+1) + "\t");
					}
					if (j == popPatToCtrl)
						Tools.printDetails(" |  ");
					Tools.printDetails(popMatrix[i][j] + "   ");
				}
				Tools.printDetails("\n");
			}
		}
		
		//check zygosity
		int e = 0;

        //changed : to = in the statement below.
		if (DETAILS) Tools.printDetails("\n\n\nZYGOSITY two by twos=", "b");
		if (DETAILS) Tools.printDetails("\n(hom) {hom dis, absent dis, hom ctr, absent ctr}\n(het) {het dis, absent dis, het ctr, absent ctr}\n(zyg) {hom dis, het dis, hom ctr, het ctr}\n");
		for (int i = 0; i < popMatrix.length; i++){
			int[] hom = {0,0}, het = {0,0}, nul = {0,0};
			int index = 0;
			for (int j = 0; j < popMatrix[0].length; j++){
				//popPatToCtrl zygosities
				if (j >= popPatToCtrl)
					index = 1;
				if (popMatrix[i][j] > 1)
					hom[index]++;
				else if ((popMatrix[i][j] == 1))
					het[index]++;
				else if (popMatrix[i][j] == 0)
					nul[index]++;
				if (j == 0)
					while (i >= hlaToNext[e+1]) // moving to next locus
						e++;
			}
			String aas="", poss="";
			boolean zyg = false;
			double pval, or;
			int[] twoby2;

			for (int m = 0; m < resArr[e][i-hlaToNext[e]].getResAAs().length; m++){
				aas += resArr[e][i-hlaToNext[e]].getResAAs()[m];
				poss += (resArr[e][i-hlaToNext[e]].getResPoss()[m]+1) + ",";
			}
			String ret = residues[e].getLocus() + "_" + aas + "-" + poss.substring(0, poss.length()-1) + ";";
			if (DETAILS) Tools.printDetails("\n" + ret, "b");
			// homozygotes associated?
			twoby2 = new int[]{hom[0], nul[0], hom[1], nul[1]};
			pval = FisherExact.exact22(twoby2) * zygosityTests;
			or = Tools.oddsRatio(twoby2);
			if (pval < Settings.GET_PVAL_CUTOFF()){
				ret += " homozygotes associated [" + Tools.formatDec(pval, Settings.GET_DEC_PLACES_PVAL()) + ", " + Tools.formatDec(or, Settings.GET_DEC_PLACES_PERC()) + "],";
			}
			else
				ret += " homozygotes NOT individually associated,";
			if (DETAILS) Tools.printDetails("\t(hom) {" + twoby2[0] + "," + twoby2[1] + "," + twoby2[2] + "," + twoby2[3] + "} p equals" + Tools.formatDec(pval, Settings.GET_DEC_PLACES_PVAL()) + ", or equals" + Tools.formatDec(or, Settings.GET_DEC_PLACES_PERC()));
			// heterozygotes associated?
			twoby2 = new int[]{het[0], nul[0], het[1], nul[1]};
			pval = FisherExact.exact22(twoby2) * zygosityTests;
			or = Tools.oddsRatio(twoby2);
			if (pval < Settings.GET_PVAL_CUTOFF()){
				ret += " heterozygotes associated [" + Tools.formatDec(pval, Settings.GET_DEC_PLACES_PVAL()) + ", " + Tools.formatDec(or, Settings.GET_DEC_PLACES_PERC()) + "],";
			}
			else
				ret += " heterozygotes NOT individually associated, ";
			if (DETAILS) Tools.printDetails("\t(het) {" + twoby2[0] + "," + twoby2[1] + "," + twoby2[2] + "," + twoby2[3] + "} p equals" + Tools.formatDec(pval, Settings.GET_DEC_PLACES_PVAL()) + ", or equals" + Tools.formatDec(or, Settings.GET_DEC_PLACES_PERC()));
			// zygosity associated?
			twoby2 = new int[]{hom[0], het[0], hom[1], het[1]};
			pval = FisherExact.exact22(twoby2) * zygosityTests;
			or = Tools.oddsRatio(twoby2);
			if (pval < Settings.GET_PVAL_CUTOFF()){
				ret += " zygosity significant [" + Tools.formatDec(pval, Settings.GET_DEC_PLACES_PVAL()) + ", " + Tools.formatDec(or, Settings.GET_DEC_PLACES_PERC()) + "]" + Settings.ARROW;
				zyg = true;
			}
			else
				ret += "no difference in zygosity.";
			if (DETAILS) Tools.printDetails("\t(zyg) {" + twoby2[0] + "," + twoby2[1] + "," + twoby2[2] + "," + twoby2[3] + "} p equals" + Tools.formatDec(pval, Settings.GET_DEC_PLACES_PVAL()) + ", or equals" + Tools.formatDec(or, Settings.GET_DEC_PLACES_PERC()));
			if (Settings.GET_PRINT_SIGNIFICANT_ONLY()){
				if (zyg) zygosities.add(ret);
			}
			else
				zygosities.add(ret);
		}
		
		if (DETAILS) Tools.printDetails("\n\n\nINTERACTIONS two by twos", "b");
		if (DETAILS) Tools.printDetails("\n(A indep B) {x1, x3, y1, y3} | {x2, x4, y2, y4}\n(B indep A) {x1, x2, y1, y2} | {x3, x4, y3, y4}\n(A inter B) {x1, x3, y1, y3}\n(A diffr B) {x2, x3, y2, y3}\n(A combd B) {x1, x4, y1, y4}\n(LD pat)    {x1, x2, x3, x4}\n(LD ctr)    {y1, y2, y3, y4}\nwhere x=dis, y=ctr, 1=++, 2=+-, 3=-+, 4=--");

        //changed : to = in the statement below.
		if (DETAILS) Tools.printDetails("\n\n----- within molecules=", "b");
		// permute matrices
		for (int m = 0; m < alleleMatrix.size(); m++){
			whichMatrix = m;
			permute(alleleMatrix.get(m));
		}
		whichMatrix = -1;

        //changed : to = in the statement below.
		if (DETAILS) Tools.printDetails("\n\n----- between molecules=", "b");
		permute(popMatrix);
	}
	
	private final void addToPopMatrix(int[][] m){
		if (m.length == 0) return;
		int x = m.length + popMatrix.length;
		int y = m[0].length;
		int[][] newMat = new int[x][y];
		for (int i = 0; i < popMatrix.length; i++){
			for (int j = 0; j < popMatrix[0].length; j++){
				newMat[i][j] = popMatrix[i][j];
			}
		}
		for (int i = popMatrix.length; i < newMat.length; i++){
			for (int j = 0; j < newMat[0].length; j++){
				newMat[i][j] = m[i-popMatrix.length][j];
			}
		}
		popMatrix = newMat;
	}
	
	public final void permute(int[][] matrix){
		int N = matrix.length;
		int k = 2;
		
		int[] s = new int[N];
		generate(s, 0, 0, k, N);
	}
	
    public final void generate(int[] s, int position, int nextInt, int k, int N) {
        if (position == k) {
        	computeCombination(s, k);
            return;
        }
        for (int i = nextInt; i < N; i++) {
            s[position] = i;
            generate(s, position + 1, i + 1, k, N);
        }
    }
	
    private void computeCombination(int[] s, int k) {
    	int col1 = s[0], col2 = s[1];
    	int hla1 = whichMatrix, hla2 = whichMatrix;
    	int add1 = 0, add2 = 0;
    	int[][] theMatrix;
    	
    	if (whichMatrix == -1){ // permuting the popMatrix
	    	for (int i = 0; i < hlaToNext.length-1; i++){
	    		if (hla1 < 0 && col1 < hlaToNext[i+1])
	    			hla1 = i;
	    		if (hla2 < 0 && col2 < hlaToNext[i+1])
	    			hla2 = i;
	    	}
	    	if (hla1 == hla2)
	    		return; // ignore combination if in popMatrix and in the same HLA (already done)
	    	
	    	theMatrix = popMatrix;
	    	add1 = hlaToNext[hla1];
	    	add2 = hlaToNext[hla2];
	    	col1 = col1 - add1;
	    	col2 = col2 - add2;
    	}
    	else{
    		theMatrix = alleleMatrix.get(whichMatrix);
    	}

        //changed == to --
    	if (DEBUG) Tools.printDetails("\n -- cols " + (col1+1+add1) + " x " + (col2+1+add2) + " of " + theMatrix.length);
    	if (resArr[hla1][col1].getResPoss()[0] != resArr[hla2][col2].getResPoss()[0] &&
    			(!ALLOW_POS_NEG_INTERACTION? resArr[hla1][col1].isPositive() == resArr[hla2][col2].isPositive() : true)){
    		int[][] twoby2info = compareColumns(col1+add1, col2+add2, theMatrix);
    		int[][] resPoss = {resArr[hla1][col1].getResPoss(), resArr[hla2][col2].getResPoss()};
    		Locus[] loci = {residues[hla1].getLocus(), residues[hla2].getLocus()};
    		char[][] aas = {resArr[hla1][col1].getResAAs(), resArr[hla2][col2].getResAAs()};
    		double[][] test = new double[(int)additionalTests][2]; // index of 0 is p-val, of 1 is OR for each of the 8 tests 
    		int x1=twoby2info[0][0],x2=twoby2info[1][0],x3=twoby2info[2][0],x4=twoby2info[3][0];
    		int y1=twoby2info[0][1],y2=twoby2info[1][1],y3=twoby2info[2][1],y4=twoby2info[3][1];
    		int[] twoby2;
    		double pval, oddsr;
    		String description = "";
    		
    		// build the identifier string for the pair
    		String posStrA = "", posStrB = "";
			String aaStrA = "", aaStrB = "";
			for (int n = 0; n < resPoss[0].length; n++){
				posStrA += (resPoss[0][n]+1) + ",";
				aaStrA += aas[0][n];
			}
			posStrA = Tools.pop(posStrA);
			for (int n = 0; n < resPoss[1].length; n++){
				posStrB += (resPoss[1][n]+1) + ",";
				aaStrB += aas[1][n];
			}
			posStrB = Tools.pop(posStrB);
			
    		String factorA = loci[0] + "_" + aaStrA + posStrA;
    		String factorB = loci[1] + "_" + aaStrB + posStrB;
    		
    		// do tests
    		LDEntry newEntry;
    		int t = 0;
    		boolean indepAB = false, indepBA = false;
    		if (DETAILS) Tools.printDetails("\n" + loci[0] + "_" + (resPoss[0][0]+1) + "-" + aas[0][0] + ", " + loci[1] + "_" + (resPoss[1][0]+1) + "-" + aas[1][0], "b");
    		
    		// A associated in B-positives?
			twoby2 = new int[]{x1, x3, y1, y3};
			pval = FisherExact.exact22(twoby2) * tests;
    		oddsr = Tools.oddsRatio(twoby2);
    		test[t][0] = (pval>1?1:pval); test[t][1] = oddsr;
    		// A associated in B-negatives?
    		twoby2 = new int[]{x2, x4, y2, y4};
    		pval = FisherExact.exact22(twoby2) * tests;
    		oddsr = Tools.oddsRatio(twoby2);
    		test[++t][0] = (pval>1?1:pval); test[t][1] = oddsr;
    		
    		// --> A independent of B
    		if (test[0][0] < Settings.GET_PVAL_CUTOFF() && test[1][0] < Settings.GET_PVAL_CUTOFF()){
    			indepAB = true;
    		}
    		if (Settings.PRINT_SIGNIFICANT_ONLY?indepAB:true){
				description += factorA + (indepAB?" ":" NOT ") + "independent of " + factorB + ", ";
    		}
    		if (DETAILS) Tools.printDetails("\n\t(A indep B) {" + x1 + "," + x3 + "," + y1 + "," + y3 + "} p equals" + Tools.formatDec(test[0][0], Settings.GET_DEC_PLACES_PVAL()) + " or equals" + Tools.formatDec(test[0][1],Settings.GET_DEC_PLACES_PERC()) + " | {" + x2 + "," + x4 + "," + y2 + "," + y4 +
                    "} p equals" + Tools.formatDec(test[1][0], Settings.GET_DEC_PLACES_PVAL()) + " or equals" + Tools.formatDec(test[1][1],Settings.GET_DEC_PLACES_PERC()));
			
    		//////////// swap residues for reverse test
    		resPoss = new int[][]{resArr[hla2][col2].getResPoss(), resArr[hla1][col1].getResPoss()};
    		loci = new Locus[]{residues[hla2].getLocus(), residues[hla1].getLocus()};
    		aas = new char[][]{resArr[hla2][col2].getResAAs(), resArr[hla1][col1].getResAAs()};
    		
    		// B associated in A-positives?
			twoby2 = new int[]{x1, x2, y1, y2};
			pval = FisherExact.exact22(twoby2) * tests;
    		oddsr = Tools.oddsRatio(twoby2);
    		test[++t][0] = (pval>1?1:pval); test[t][1] = oddsr;
    		// B associated in A-negatives?
    		twoby2 = new int[]{x3, x4, y3, y4};
    		pval = FisherExact.exact22(twoby2) * tests;
    		oddsr = Tools.oddsRatio(twoby2);
    		test[++t][0] = (pval>1?1:pval); test[t][1] = oddsr;
    		
    		// --> B independent of A
    		if (test[2][0] < Settings.GET_PVAL_CUTOFF() && test[3][0] < Settings.GET_PVAL_CUTOFF()){
    			indepBA = true;
    		}
    		if (Settings.PRINT_SIGNIFICANT_ONLY?indepBA:true){
				description += factorB + (indepBA?" ":" NOT ") + "independent of " + factorA + ", ";
    		}
    		if (DETAILS) Tools.printDetails("\n\t(B indep A) {" + x1 + "," + x2 + "," + y1 + "," + y2 + "} p equals" + Tools.formatDec(test[2][0], Settings.GET_DEC_PLACES_PVAL()) + " or equals" + Tools.formatDec(test[2][1],Settings.GET_DEC_PLACES_PERC()) + " | {" + x3 + "," + x4 + "," + y3 + "," + y4 +
                    "} p equals" + Tools.formatDec(test[3][0], Settings.GET_DEC_PLACES_PVAL()) + " or equals" + Tools.formatDec(test[3][1],Settings.GET_DEC_PLACES_PERC()));
    		
    		boolean inter = false;
    		// --> A and B show interaction
    		if (test[0][0] < Settings.GET_PVAL_CUTOFF() && test[2][0] < Settings.GET_PVAL_CUTOFF()){
    			inter = true;
    		}
    		if (Settings.PRINT_SIGNIFICANT_ONLY?inter:true){
    			description += (inter?"":"DO NOT ") + "interact, ";
    		}
    		if (DETAILS) Tools.printDetails("\n\t(A inter B) {" + x1 + "," + x3 + "," + y1 + "," + y3 + "} p equals" + Tools.formatDec(test[0][0], Settings.GET_DEC_PLACES_PVAL()) + " or equals" + Tools.formatDec(test[0][1],Settings.GET_DEC_PLACES_PERC()) + " | {" + x1 + "," + x2 + "," + y1 + "," + y2 +
                    "} p equals" + Tools.formatDec(test[2][0], Settings.GET_DEC_PLACES_PVAL()) + " or equals" + Tools.formatDec(test[2][1],Settings.GET_DEC_PLACES_PERC()));
    		
    		// to make sure the residues are named
    		if (Settings.PRINT_SIGNIFICANT_ONLY && !indepBA && !indepAB)
    			description += factorA + ", " + factorB + " ";
    		
    		// Difference between A and B associations?
    		twoby2 = new int[]{x2, x3, y2, y3};
    		pval = FisherExact.exact22(twoby2) * tests;
    		oddsr = Tools.oddsRatio(twoby2);
    		test[++t][0] = (pval>1?1:pval); test[t][1] = oddsr;
    		
    		boolean diff = false;
    		// --> Difference btw A and B
    		if (test[t][0] < Settings.GET_PVAL_CUTOFF()){
    			diff = true;
    		}
    		if (Settings.PRINT_SIGNIFICANT_ONLY?diff:true){
    			description += "associations" + (diff?"":" DO NOT") + " differ, ";
    		}
    		if (DETAILS) Tools.printDetails("\n\t(A diffr B) {" + x2 + "," + x3 + "," + y2 + "," + y3 + "} p equals" + Tools.formatDec(test[4][0], Settings.GET_DEC_PLACES_PVAL()) + " or equals" + Tools.formatDec(test[4][1],Settings.GET_DEC_PLACES_PERC()));
		
    		
    		// ADDITIONAL TESTS
    		// --> Combined A-B association
    		twoby2 = new int[]{x1, x4, y1, y4};
    		pval = FisherExact.exact22(twoby2) * additionalTests;
    		oddsr = Tools.oddsRatio(twoby2);
    		test[++t][0] = (pval>1?1:pval); test[t][1] = oddsr;
    		boolean comb = false;
    		if (test[t][0] < Settings.GET_PVAL_CUTOFF()){
    			comb = true;
    		}
    		if (Settings.PRINT_SIGNIFICANT_ONLY?comb:true){
    			description += (comb?"":"DO NOT ") + "have combined action, ";
    		}
    		if (DETAILS) Tools.printDetails("\n\t(A combd B) {" + x1 + "," + x4 + "," + y1 + "," + y4 + "} p equals" + Tools.formatDec(test[5][0], Settings.GET_DEC_PLACES_PVAL()) + " or equals" + Tools.formatDec(test[5][1],Settings.GET_DEC_PLACES_PERC()));
    		// --> LD in pat
    		twoby2 = new int[]{x1, x2, x3, x4};
    		pval = FisherExact.exact22(twoby2) * additionalTests;
    		oddsr = Tools.oddsRatio(twoby2);
    		test[++t][0] = (pval>1?1:pval); test[t][1] = oddsr;
    		boolean ldp = false;
    		if (test[t][0] < Settings.GET_PVAL_CUTOFF()){
    			ldp = true;
    		}
    		if (Settings.PRINT_SIGNIFICANT_ONLY?ldp:true){
    			description += (ldp?"":"NOT ") + "in LD (CASE), ";
    		}
    		if (DETAILS) Tools.printDetails("\n\t(LD pat)    {" + x1 + "," + x2 + "," + x3 + "," + x4 + "} p equals" + Tools.formatDec(test[6][0], Settings.GET_DEC_PLACES_PVAL()) + " or equals" + Tools.formatDec(test[6][1],Settings.GET_DEC_PLACES_PERC()));
    		// --> LD in ctrl
    		twoby2 = new int[]{y1, y2, y3, y4};
    		pval = FisherExact.exact22(twoby2) * additionalTests;
    		oddsr = Tools.oddsRatio(twoby2);
    		test[++t][0] = (pval>1?1:pval); test[t][1] = oddsr;
    		boolean ldc = false;
    		if (test[t][0] < Settings.GET_PVAL_CUTOFF()){
    			ldc = true;
    		}
    		if (Settings.PRINT_SIGNIFICANT_ONLY?ldc:true){
    			description += (ldc?"":"NOT ") + "in LD (CTRL), ";
    		}
    		if (DETAILS) Tools.printDetails("\n\t(LD ctr)    {" + y1 + "," + y2 + "," + y3 + "," + y4 + "} p equals" + Tools.formatDec(test[7][0], Settings.GET_DEC_PLACES_PVAL()) + " or equals" + Tools.formatDec(test[7][1],Settings.GET_DEC_PLACES_PERC()));

    		newEntry = new LDEntry(loci, resPoss, aas, test, description);
        	patCtrlEntries.add(newEntry);
    	}
    }
    
    private int[][] compareColumns(int col1, int col2, int[][] matrix){
    	int[] c1 = {0,0}, c2 = {0,0}, both = {0,0}, none = {0,0}, missing = {0,0};
    	int index = 0; // 0-patients 1-controls
    	int patCtrl;
    	if (whichMatrix == -1) // permuting the popMatrix
    		patCtrl = popPatToCtrl;
    	else
    		patCtrl = patToCtrl[whichMatrix];
    	
    	for (int i = 0; i < matrix[col1].length; i++){
    		if (i >= patCtrl)
    			index = 1;
    		
    		// get the pop count
    		int popCount = 1;
    		if (whichMatrix > -1){
    			if (index == 0)
    				popCount = (pats[whichMatrix].getAllelesMap().get(pats[whichMatrix].getAlleles()[i])[2]).intValue();
    			else
    				popCount = (ctrls[whichMatrix].getAllelesMap().get(ctrls[whichMatrix].getAlleles()[i-patToCtrl[whichMatrix]])[2]).intValue();
    		}
    		
    		if      (matrix[col1][i] > 0 && matrix[col2][i] > 0)
    			both[index] += popCount;
    		else if (matrix[col1][i] > 0 && matrix[col2][i] == 0)
    			  c1[index] += popCount;
    		else if (matrix[col1][i] == 0 && matrix[col2][i] > 0)
    			  c2[index] += popCount;
    		else if (matrix[col1][i] == 0 && matrix[col2][i] == 0)
    			none[index] += popCount;
    		else
    			missing[index]++;
    	}
    	if (DEBUG){

            //changed : to = in the statement below.
    		Tools.printDetails("\n" + " pat (++)x1=" + both[0] + " (+-)x2=" + c1[0] + " (-+)x3=" + c2[0] + " (--)x4=" + none[0] + (missing[0]>0?" -missing="+missing[0]:""));
    		Tools.printDetails("\n" + " ctr (++)y1=" + both[1] + " (+-)y2=" + c1[1] + " (-+)y3=" + c2[1] + " (--)y4=" + none[1] + (missing[1]>0?" -missing="+missing[1]:""));
    	}
    	return new int[][]{both, c1, c2, none};
    }
    
    public class LDEntry{
    	protected int[][] positions;
    	protected char[][] aas;
    	protected double[][] tests;
    	protected String description;
    	protected Locus[] hla;
    	
    	public LDEntry(Locus[] hla, int[][] positions, char[][] aas, double[][] tests, String description){
    		this.hla = hla;
    		this.positions = positions;
    		this.aas = aas;
    		this.tests = tests;
    		this.description = description;
    	}

    	public final Locus[] getHla(){
    		return hla;
    	}
    	
		public final char[][] getAAs(){
			return aas;
		}

		public final double getTestPval(int t){
			return tests[t][0];
		}
		
		public final double getTestOR(int t){
			return tests[t][1];
		}
		
		public final int getTestCount(){
			return tests.length;
		}

		public final int[][] getPositions(){
			return positions;
		}
		
		public final String getDescription(){
			return description;
		}
    }
    
    public String toString(){
    	StringBuffer sb = new StringBuffer();
    	if (resArr.length > 0){
    		
    		// zygosity
	    	Iterator<String> zygit = zygosities.iterator();

            //changed : to = in the statement below.
	    	sb.append("\n\n" + Tools.underline("Zygosity analysis for " + pats[0].getTitle() + " and " + ctrls[0].getTitle() + "="));
	    	sb.append("\n> p-value correction is " + (int)zygosityTests + " (equal to zygosity tests)");
	    	if (zygosities.size() > 0){
	    		sb.append("\nDescription [p-val, OR]");
	    	}
	    	while (zygit.hasNext()){
	    		sb.append("\n" + zygit.next());
	    	}
	    	if (zygosities.size() == 0)
	    		sb.append("\nNo significance in zygosity found.");
	    	
	    	// interactions
	    	if (patCtrlEntries.size() > 1){ // if we have 2 or more significant residues
	    		// sort
	    		TreeMap<String, LDEntry> tm = new TreeMap<String, LDEntry>();
		    	Iterator<LDEntry> it = patCtrlEntries.iterator();
		    	int id = 0;
		    	while (it.hasNext()){
		    		LDEntry curr = it.next();
		    		String poss = String.valueOf(curr.getPositions()[0][0]);
		    		for (int p = String.valueOf(curr.getPositions()[0][0]).length(); p < 3; p++)
		    			poss = "0" + poss;
		    		String str_ld = curr.getHla()[0] + poss + curr.getAAs()[0][0] + "" + id++;
		    		tm.put(str_ld, curr);
		    	}
		    	Iterator<LDEntry> tmit = tm.values().iterator();


                //changed : to = in the statement below.
		    	sb.append("\n\n" + Tools.underline("=Interaction analysis for " + pats[0].getTitle() + " and " + ctrls[0].getTitle() + "="));

                //changed = to equal to
		    	sb.append("\n> p-value correction is " + (int)tests + " (equal to tests 1-5, for strongest association), " + (int)additionalTests + " (for less critical tests, 6-8)");
		    	int count = 0;
		    	while (tmit.hasNext()){
		    		LDEntry curr = tmit.next();
	    			if (count == 0)
	    				sb.append("\nDescription [test; p-val, OR]");
					sb.append("\n" + Tools.pop(curr.getDescription(), 2) + ".\n\t");
					for (int t = 0; t < curr.getTestCount(); t++){
						if (Settings.PRINT_SIGNIFICANT_ONLY?curr.getTestPval(t) < Settings.GET_PVAL_CUTOFF():true)
							sb.append("[" + (t+1) + "; " + Tools.formatDec(curr.getTestPval(t),Settings.GET_DEC_PLACES_PVAL()) + ", " + Tools.formatDec(curr.getTestOR(t),Settings.GET_DEC_PLACES_PERC()) + "]" + (curr.getTestPval(t)<Settings.GET_PVAL_CUTOFF() && curr.getTestOR(t)>0? Settings.ARROW:"") + "\t");
					}
		    		count++;
		    	}
		    	if (count == 0) sb.append("\nNo interactions found.");
		    }
    	}
	    sb.append(Settings.EMPTY_STRING + Settings.EMPTY_STRING);
	    return sb.toString();
    }
}
