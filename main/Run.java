package main;

import utils.Settings;
import utils.Tools;
import web.HLAAlign;

import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.Calendar;
import java.util.LinkedList;

public class Run {
	public static PrintStream LOG;
	
	public static void main(String[] args){
		try{
			long start = System.currentTimeMillis();
			Settings.load();
			
			FileOutputStream foutLog = new FileOutputStream(Settings.WORKING_DIR + "log.txt");
			LOG = new PrintStream(foutLog);
			Tools.printAll("SKDM HLA Tools " + Settings.VERSION + " by S.Kanterakis & D.Monos, UPENN/CHOP (c)2007" + "\nTimestamp: " + Calendar.getInstance().getTime().toString() + "\n\n");
			
			ParseFile pat = null, ctrl = null;
			if (args == null){
				File f1 = new File(Settings.WORKING_DIR + "pat.txt");
				File f2 = new File(Settings.WORKING_DIR + "ctrl.txt");
				if (f1.exists() && f1.length() > 0){
					pat = new ParseFile(f1, "CASES");
					if (f2.exists() && f2.length() > 0){
						ctrl = new ParseFile(f2, "CONTROLS");
					}
				}
				else if (f2.exists() && f2.length() > 0){
					pat = new ParseFile(f2, "CONTROLS");
				}
				else{
					Tools.print("No input found!", "font", "color=red");
					Tools.print("\nPlease paste a case and/or a control dataset on the left,\nor type something like:\nA\n0101");
					return;
				}
			}
			else if (args.length == 2){
				if (args[0] != null && args[1] != null){
					pat = new ParseFile(args[0], "CASES");
					ctrl = new ParseFile(args[1], "CONTROLS");
					FileOutputStream fout;
				    fout = new FileOutputStream(Settings.WORKING_DIR + "pat.txt");
				    new PrintStream(fout).println(args[0]);
				    fout.close();
				    fout = new FileOutputStream(Settings.WORKING_DIR + "ctrl.txt");
				    new PrintStream(fout).println(args[1]);
				    fout.close();
				}
				else if (args[0] == null){
					pat = new ParseFile(args[1], "CONTROLS");
					FileOutputStream fout;
				    fout = new FileOutputStream(Settings.WORKING_DIR + "pat.txt");
				    fout.close();
				    fout = new FileOutputStream(Settings.WORKING_DIR + "ctrl.txt");
				    new PrintStream(fout).println(args[1]);
				    fout.close();
				}
				else if (args[1] == null){
					pat = new ParseFile(args[0], "CASES");
					FileOutputStream fout;
				    fout = new FileOutputStream(Settings.WORKING_DIR + "pat.txt");
				    new PrintStream(fout).println(args[0]);
				    fout.close();
				    fout = new FileOutputStream(Settings.WORKING_DIR + "ctrl.txt");
				    fout.close();
				}
				else
					return;
				
			}

			
			LinkedList<Delta> del = new LinkedList<Delta>();
			ListUnique[] lu_ctrl = null, lu_pat = null;
			
			if (ctrl != null){
				Tools.printAll("\n" + pat.getTitle() + " vs " + ctrl.getTitle(), "span", "style='background-color:black; color:white; font-weight:bold;'");
			}


            //changed : to = in the statement below.
			if (pat.getHLAs() == null)
				return;
			Tools.print(pat.titleToString(), "b");
			lu_pat = new ListUnique[pat.getHLAs().size()];
			for (int i = 0; i <pat.getHLAs().size(); i++){
				lu_pat[i] = new ListUnique(pat.getHLAs().get(i), pat.getAlleles().get(i), pat.getPopCounts()[i], pat.getTitle());
				String strLu = lu_pat[i].toString();
				Tools.print(strLu.split("=")[0],"b");
				Tools.print(strLu.split("=")[1]);

			}

            //changed : to = in the statement below.
			if (ctrl != null){
				Tools.print(ctrl.titleToString(), "b");
				lu_ctrl = new ListUnique[ctrl.getHLAs().size()];
				for (int i = 0; i <ctrl.getHLAs().size(); i++){
					lu_ctrl[i] = new ListUnique(ctrl.getHLAs().get(i), ctrl.getAlleles().get(i), ctrl.getPopCounts()[i], ctrl.getTitle());
					String strLu = lu_ctrl[i].toString();
                    Tools.print(strLu.split("=")[0], "b");
                    Tools.print(strLu.split("=")[1]);
				}
				for (int i = 0; i <pat.getHLAs().size(); i++){
					if (ctrl.getHLAs().contains(pat.getHLAs().get(i))){
						del.add(new Delta(lu_ctrl[ctrl.getHLAs().indexOf(pat.getHLAs().get(i))], lu_pat[i], ctrl.getHLAs().get(i)));
						String strDel = del.get(i).toString();
						Tools.print(strDel.split("=")[0], "b");
						Tools.print(strDel.split("=")[1]);

					}else
						if (!Settings.SUPPRESS_WARNINGS)Tools.print("\n\n*Locus HLA-" + pat.getHLAs().get(i) + " is not contained in both " + pat.getTitle() + " and " + ctrl.getTitle() + "\nand is therefore not eligible for Delta analysis!", "font", "color=red");
				}
				for (int i = 0; i <ctrl.getHLAs().size(); i++){
					if (!pat.getHLAs().contains(ctrl.getHLAs().get(i)))
						if (!Settings.SUPPRESS_WARNINGS)Tools.print("\n\n*Locus HLA-" + ctrl.getHLAs().get(i) + " is not contained in both " + pat.getTitle() + " and " + ctrl.getTitle() + "\nand is therefore not eligible for Delta analysis!", "font", "color=red");
				}
			}
			if (del.size() > 0){
				HLAAlign[] aligns = new HLAAlign[del.size()];
				Residue[] res = new Residue[aligns.length];
                //changed : to = in the statement below.
				Tools.printDetails("\n\nRESIDUE two by twos=","b");
				Tools.printDetails("\n{patient+, patient-, contol+, contol-}");
				for (int i = 0; i <del.size(); i++){
					aligns[i] = new HLAAlign(del.get(i).loc, del.get(i).getAllelesSorted(), del.get(i).getDeltas(), del.get(i).getPvals(), del.get(i).getORs());
					String strAl = aligns[i].toString();
                    Tools.print(strAl.split("=")[0], "b");
                    Tools.print(strAl.split("=")[1]);

					res[i] = new Residue(aligns[i], ctrl.getAllelesOfLoc(ctrl.getHLAs().get(i)), ctrl.getPopCounts()[i], pat.getAllelesOfLoc(pat.getHLAs().get(i)), pat.getPopCounts()[i]);
					String strRes = res[i].toString();
                    //tring[] strResArr = strRes.split("=");
                    Tools.print(strRes.split("=")[0], "b");
                    Tools.print(strRes.split("=")[1]);
                    Tools.print(strRes.split("=")[2], "b");
                    Tools.print(strRes.split("=")[3]);
				}
				LinkageDisequilibrium ld = new LinkageDisequilibrium(lu_pat,lu_ctrl,res);
				String strLD = ld.toString();

                //String[] strLDArr = strLD.split("=");
                Tools.print(strLD.split("=")[0], "b");
                Tools.print(strLD.split("=")[1]);
                Tools.print(strLD.split("=")[2], "b");
                Tools.print(strLD.split("=")[3]);
			}
			else{
				int pn = pat.getHLAs().size();
				int cn = ctrl != null? ctrl.getHLAs().size() : 0;
				HLAAlign[] aligns = new HLAAlign[cn + pn];
				for (int i = 0; i < pn; i++){

                    //String[] strLDArr = strLD.split("=");
					aligns[i] = new HLAAlign(pat.getHLAs().get(i), lu_pat[i].getAllelesList());
					String strAl = aligns[i].toString();
                    Tools.print(strAl.split("=")[0], "b");
                    Tools.print(strAl.split("=")[1]);
				}
				for (int i = pn; i < pn+cn; i++){

                    //String[] strLDArr = strLD.split("=");
					aligns[i] = new HLAAlign(ctrl.getHLAs().get(i-pn), lu_ctrl[i-pn].getAllelesList());
					String strAl = aligns[i].toString();
                    Tools.print(strAl.split("=")[0], "b");
                    Tools.print(strAl.split("=")[1]);
				}
			}
			long finish = System.currentTimeMillis() - start;
			Tools.print("\n\n\nAnalysis completed in " + Tools.formatDec((int)finish, 3) + " seconds." + "\n*EOF*****************************************************\n");
			foutLog.close();
		}
		catch(Exception e){
			Tools.printErr("\n\n!!!\nThere was an error during execution: " + e.toString() + "\n!!!\n\n");
			e.printStackTrace();
		}
	}
}
