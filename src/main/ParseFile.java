package src.main;

import src.utils.Settings;
import src.utils.Tools;

import java.io.*;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.StringTokenizer;

public class ParseFile {
	private LinkedList<Locus> HLAs;
	private int[] popCounts;
	private int loci;
	private LinkedList<String> IDs;
	private LinkedList<LinkedList<String>> alleles;
	private String title;
	private boolean noIDs = true;
	
	public ParseFile(File f, String defaultTitle){
		HLAs = null;
		BufferedReader br;
		try{
			br = new BufferedReader(new FileReader(f));
		}catch (FileNotFoundException e){
			Tools.print("\nCould not read input file. Aborting...", "font", "color=red");
			return;
		}
		continueParse(br, f.getAbsolutePath(), defaultTitle);
	}
	
	public ParseFile(String input, String defaultTitle){
		BufferedReader br;
        try{
			br = new BufferedReader(new StringReader(input));
		}catch (Exception e){
			Tools.print("\nCould not read input data. Aborting...", "font", "color=red");
			return;
		}
		continueParse(br, "(textbox input)", defaultTitle);
	}
	
	private void continueParse (BufferedReader br, String inTitle, String defaultTitle){
		String line, line2;
		StringTokenizer st;
		try {
			line = br.readLine();
			while (line != null && (line.length()<1 || line.matches("[ \t]+"))) line = br.readLine();
			if (line == null){
				Tools.print("\nInput file is empty. Aborting...", "font", "color=red");
				return;
			}
			title = line;
			line = br.readLine(); line2 = line;
			boolean noTitle = false;
			
			// determine if title exists
			if (line != null && Character.isDigit(line.charAt(0))){
				line = title;
				noTitle = true;
				title = defaultTitle;
			}
			else{
				title = title.replaceAll("\t"," ");
			}
			
			// try to get allele headers
			while (line != null && (line.length()<1 || line.matches("[ \t]+"))) line = br.readLine();
			if (line == null){
				Tools.print("\nInput does not contain any HLA information. Aborting...", "font", "color=red");
				Tools.print("\nPlease paste a case and/or a control dataset on the left,\nor type something like:\nA\n0101");
				return;
			}
			
			st = new StringTokenizer(line, "\t", true);
			String loc = "", tok = "";
			int toks = 0; int locERR = 0;
			HLAs = new LinkedList<Locus>();
			if (st.hasMoreTokens()) tok = st.nextToken();
			boolean theEnd = !st.hasMoreTokens();
			do{
				if (!noIDs && toks == 1) if (st.hasMoreTokens()) tok = st.nextToken();
				if (tok.equals("\t")){
					if (toks == 0){ // in case id header is missing
						noIDs = false; 
						toks++;
					}
					else{
						if (st.hasMoreTokens()) tok = st.nextToken();
						if (!tok.equals("\t")){ // user has a header on each column
							if (st.hasMoreTokens()) tok = st.nextToken();
							else
								theEnd = true;
						}
						if (st.hasMoreTokens()) tok = st.nextToken();
						else
							theEnd = true;
					}
				}
				else{
					loc = tok.toUpperCase().replaceAll("HLA","").replaceAll("[^A-Z0-9]","");
					
					/*if (loc.startsWith(Locus.A.toString()))
						HLAs.add(Locus.A);
					else if (loc.startsWith("B"))
						HLAs.add(Locus.B);
					else if (loc.startsWith("Cw") || loc.startsWith("C"))
						HLAs.add(Locus.Cw);
					else if (loc.startsWith("DRB1") || loc.startsWith("DRB"))
						HLAs.add(Locus.DRB1);
					else if (loc.startsWith("DQB1") || loc.startsWith("DQB"))
						HLAs.add(Locus.DQB1);
					else if (loc.startsWith("DPB1") || loc.startsWith("DPB"))
						HLAs.add(Locus.DPB1);
					else if (loc.startsWith("DQA1") || loc.startsWith("DQA"))
						HLAs.add(Locus.DQA1);
					else if (loc.startsWith("DPA1") || loc.startsWith("DPA"))
						HLAs.add(Locus.DPA1);
					else if (loc.startsWith("MICA"))
						HLAs.add(Locus.MICA);
					else if (loc.startsWith("MICB"))
						HLAs.add(Locus.MICB);
                    else if (loc.equals("") || loc.matches("[ ]+"));
					else if (toks == 0)
						noIDs = false; // file contains identifiers*/
                    if (loc.startsWith(Locus.A.toString()))
                        HLAs.add(Locus.A);
                    else if (loc.startsWith(Locus.B.toString()))
                        HLAs.add(Locus.B);
                    else if (loc.startsWith(Locus.Cw.toString()))
                        HLAs.add(Locus.Cw);
                    else if (loc.startsWith(Locus.C.toString()))
                        HLAs.add(Locus.C);
                    else if (loc.startsWith("DRB"))
                        HLAs.add(Locus.DRB1);
                    else if (loc.startsWith("DQB"))
                        HLAs.add(Locus.DQB1);
                    else if (loc.startsWith("DPB"))
                        HLAs.add(Locus.DPB1);
                    else if (loc.startsWith("DQA"))
                        HLAs.add(Locus.DQA1);
                    else if (loc.startsWith("DPA"))
                        HLAs.add(Locus.DPA1);
                    else if (loc.startsWith("MICA"))
                        HLAs.add(Locus.MICA);
                    else if (loc.startsWith("MICB"))
                        HLAs.add(Locus.MICB);
                    else if (loc.startsWith("DMA"))
                        HLAs.add(Locus.DMA1);
                    else if (loc.startsWith("DMB"))
                        HLAs.add(Locus.DMB1);
                    else if (loc.startsWith("DOA"))
                        HLAs.add(Locus.DOA);
                    else if (loc.startsWith("DOB"))
                        HLAs.add(Locus.DOB);
                    else if (loc.startsWith("DPB"))
                        HLAs.add(Locus.DPB1);
                    else if (loc.startsWith("DQB"))
                        HLAs.add(Locus.DQB1);
                    else if (loc.startsWith("DRA"))
                        HLAs.add(Locus.DRA);
                    else if (loc.startsWith("E"))
                        HLAs.add(Locus.E);
                    else if (loc.startsWith("F"))
                        HLAs.add(Locus.F);
                    else if (loc.startsWith("G"))
                        HLAs.add(Locus.G);
                    else if (loc.startsWith(Locus.TAP1.toString()))
                        HLAs.add(Locus.TAP1);
                    else if (loc.startsWith(Locus.TAP2.toString()))
                        HLAs.add(Locus.TAP2);
                    else if (loc.equals("") || loc.matches("[ ]+"));
                    else if (toks == 0)
                        noIDs = false; // file contains identifiers
                    else{
						HLAs.add(Locus.ERR);
						if (!Settings.SUPPRESS_WARNINGS && toks>0) Tools.print("\nLocus '" + loc + "' cannot be recognized and will be ignored!\n*Please use full locus names (ex. DRB1 instead of DR)", "font", "color=red");
						locERR++;
					}
					toks++;
					if (st.hasMoreTokens()) tok = st.nextToken();
					else
						theEnd = true;
				}
			}while(!theEnd);
			loci = HLAs.size();
			if (loci == locERR){
				Tools.print("\n**No valid HLA locus headers found in input.\n  Please annotate each set of columns with an HLA locus header (ex. HLA-A)\n", "font", "color=red");
				//return;
			}
			
			// initialize allele arrays
			alleles = new LinkedList<LinkedList<String>>();
			for (int i = 0; i < loci; i++){
				alleles.add(new LinkedList<String>());
			}
			
			// adjust for no title
			if (noTitle)
				line = line2;
			else
				line = br.readLine();
			
			IDs = new LinkedList<String>();
			int row = 0, col = 0;
			boolean isTab = false, flag = false;
			String token = null;
			while (line != null && (line.length()<1 || line.matches("[ \t]+"))) line = br.readLine();
			while (line != null){				
				st = new StringTokenizer(line, "\t", true);
				if (st.countTokens() > loci*4+1 || st.countTokens() < (noIDs? 1 : 2)){
					if (!Settings.SUPPRESS_WARNINGS)Tools.print("\nProblem with line \"" + line + "\". Moving on to next line...\n", "font", "color=blue");
				}
				else{
					if (noIDs) col = 1;
					else col = 0;
					token = null; flag = false; isTab = false;
					while (st.hasMoreTokens() || flag){
						if (!flag) {
							token = st.nextToken();
							if (noIDs && (!token.equals("\t") && (token.length() < 2 || !Character.isDigit(token.charAt(0))))){
								if (!Settings.SUPPRESS_WARNINGS)Tools.print("\nWarning: allele '" + token + "' is not a valid HLA allele and was removed!", "font", "color=blue");
								token = Settings.BLANK;
							}
						}
						else if (flag) {
							token = "\t";
							flag = false;
						}
						
						if (!isTab && token.equals("\t") && col > 0){
							flag = true;
							token = Settings.BLANK;
						}
						
						if (token.equals("\t") && col > 0){
							//wait for next token
						}
						else if (col == 0){
							if (token.equals("\t")){
								Tools.print("\n*Row identifier missing from row " + (row+1) +"...", "font", "color=blue");
								IDs.add("*noid*");
								flag = true;
							}else
							IDs.add(token);
						}
						else if (col % 2 == 1){
							alleles.get((col+1)/2-1).add(token);
						}
						else{
							alleles.get(col/2-1).add(token);
						}
						if (!isTab) col++;
						isTab = !isTab;
						
						if (!st.hasMoreTokens() && alleles.get(col/2-1).size() < (row+1)*2){
							alleles.get(col/2-1).add(Settings.BLANK);
						}
					}
					row++;
				}
				line = br.readLine();
				while (line != null && (line.length()<1 || line.matches("[ \t]+"))) line = br.readLine();
			}
		} catch (Exception e) {
			Tools.print("\nAn error occured: " + e.getCause(), "font", "color=red");
			e.printStackTrace();
			return;
		}
		popCounts = new int[alleles.size()];
		clearHomozygotes();
	}
	
	private void clearHomozygotes(){
		String s1,s2;
		for (int hla = 0; hla < alleles.size(); hla++){
			for (int i = 0; i < alleles.get(hla).size(); i=i+2){
				s1 = alleles.get(hla).get(i); s2 = alleles.get(hla).get(i+1);
				if (s1.equals(Settings.BLANK)) {alleles.get(hla).set(i,""); s1="";}
				if (s2.equals(Settings.BLANK)) {alleles.get(hla).set(i+1,""); s2="";}
				if (!s1.equals("") || !s2.equals("")){
					if (s1.equals(s2) || (s1.length()>0 && s2.startsWith(s1)))
						alleles.get(hla).set(i+1,"");
					else if (s2.length()>0 && s1.startsWith(s2)){
						alleles.get(hla).set(i,s2);
						alleles.get(hla).set(i+1,"");
					}
					popCounts[hla]++;
				}
			}
			//System.out.println(HLAs.get(hla) + ": " + popCounts[hla]);
		}
	}
	
	public ArrayList<String> getAllelesOfLoc(Locus X){
		ArrayList<String> ret = null;
		if (HLAs.contains(X)){
			int index = HLAs.indexOf(X);
			ret = new ArrayList<String>(alleles.get(0).size());
			for (int i = 0; i < alleles.get(index).size(); i++)
				ret.add(i, alleles.get(index).get(i));
		}
		return ret;
	}
	
	public LinkedList<Locus> getHLAs(){
		return HLAs;
	}
	
	public LinkedList<LinkedList<String>> getAlleles(){
		return alleles;
	}
	
	public int[] getPopCounts(){
		return popCounts;
	}
	
	public String titleToString(){
		return "\n\n\n" + Tools.underline(title);
	}
	
	public String getTitle(){
		return title;
	}
	
	public String toString(){
		StringBuffer sb = new StringBuffer();
		sb.append("\n\n" + Tools.underline(title) + "\n");
		if (!noIDs)
			sb.append("\t");
		for (int i = 0; i < loci; i++){
			sb.append("HLA-" + HLAs.get(i) + ":");
			if (i < loci-1)
				sb.append("\t\t");
		}
		if (alleles != null && alleles.size() > 0){
			for (int k = 0; k < alleles.get(0).size(); k=k+2){
				sb.append("\n" + (noIDs? "" : (IDs.get(k/2) + "\t")));
				for (int i = 0; i < loci; i++){
					sb.append(alleles.get(i).get(k)+ "\t" + alleles.get(i).get(k+1)+ (i==loci-1?"":"\t"));
				}
			}
			sb.append("\n" + (alleles.get(0).size()/2) + "\t" + (noIDs? "lines":"identifiers") + " present.");
			return sb.toString();
		}
		return "";
	}
}