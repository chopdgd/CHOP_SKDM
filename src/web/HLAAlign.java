package src.web;

import edu.chop.dgd.process.FileDownloader;
import src.main.Locus;
import org.apache.commons.lang3.EnumUtils;
import src.utils.Settings;
import src.utils.Tools;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.net.HttpURLConnection;
import java.net.URL;
import java.net.URLConnection;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

@SuppressWarnings("unchecked")
public class HLAAlign {
	private LinkedList<String> allelesAsc;
	private ArrayList<String> allelesRetrieved;
	private ArrayList<String> allelesIn;
    private ArrayList<String> allelesTobeRemoved;
	private ArrayList<String> deltas = null;
	private ArrayList<Double> pvals = null;
	private ArrayList<Double> ors = null;
	private Locus loc;
	private String ref;
	private char[][] aligns;
	private int allele_count;
	private char[] refseq;

	
	private static final int COLS = 10;
	private static final int COL = 100;
	private static final int MAX_HEADER = 10;
	private static final boolean DEBUG = false;
	private static final boolean ONLINE = false;
    private static final boolean DOWNLOAD = true;

    //public static final int ALIGNMENTS_PRINT_UNTIL_POS = 155;
	public static int ALIGNMENTS_PRINT_UNTIL_POS = 300;
	
	public ArrayList<String> getAlleles(){
		return allelesIn;
	}
	
	public Locus getLocus(){
		return loc;
	}
	
	public String getReference(){
		return ref;
	}
	
	public char[][] getAligns(){
		return aligns;
	}
	
	public ArrayList<Double> getPvals(){
		return pvals;
	}
	
	public ArrayList<Double> getORs(){
		return ors;
	}
	
	public int getNegativeBreak(){
		boolean found = false;
		int i = 0;
		int negi = -1;
		while (!found && i < deltas.size()){
			if (deltas.get(i).charAt(0) == '-'){
				negi = i;
				found = true;
			}
			i++;
		}
		return negi;
	}
	
	public char[] getRefSeq(){
		return refseq;
	}
	
	public int getRefSeqIndex(){
		boolean found = false;
		int i = 0;
		int refi = -1;
		while (!found && i < allelesIn.size()){
			if (ref.startsWith(allelesIn.get(i))){
				refi = i;
				found = true;
			}
			i++;
		}
		return refi;
	}
	
	public HLAAlign(Locus loc, ArrayList<String> allelesIn) throws Exception {
		//this(loc, new ArrayList(Arrays.asList(allelesInArr)), null, null);
		this(loc, allelesIn, null, null, null);
	}

	public HLAAlign(Locus loc, ArrayList<String> allelesIn, ArrayList<String> deltas, ArrayList<Double> pvals, ArrayList<Double> ors) throws Exception {
		this.allelesIn = allelesIn;
		if (deltas != null && deltas.size() == allelesIn.size()){
			this.deltas = deltas;
			if (pvals != null && pvals.size() == allelesIn.size()) this.pvals = pvals;
			if (ors != null && ors.size() == allelesIn.size()) this.ors = ors;
		}
		allelesAsc = new LinkedList();
		TreeSet <String>ts = new TreeSet();
		for (int i = 0; i<allelesIn.size(); i++){
			ts.add(allelesIn.get(i));
		}
		String first;
		for (int i = 0; i<allelesIn.size(); i++){    
			first = ts.first();
			allelesAsc.add(first);
			ts.remove(first);
		}
		this.loc = loc;
		ref = loc.getRefSeq();
		
		allele_count = allelesAsc.size();
		allelesRetrieved = new ArrayList<String>(allele_count);
		aligns = new char[allele_count][COL * COLS];
		refseq = new char[COL * COLS];
		getAlign();
		sortAligns();
	}
	
	public void getAlign() throws Exception {
		String response = null; String modResponse = null;

        //DEBUG doesn't work..
        if (DEBUG)
            response = "<!DOCTYPE html PUBLIC '-//W3C//DTD XHTML 1.0 Transitional//EN' 'http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd'><html xmlns='http://www.w3.org/1999/xhtml' lang='eng'><!-- InstanceBegin template='/Templates/new_template_no_menus_max.dwt' codeOutsideHTMLIsLocked='false' --><head><meta http-equiv='Content-Type' content='text/html; charset=iso-8859-1' /><meta name='Keywords' content='bioinformatics, IMGT, immunogenetics, database, HLA, submission of sequences, ebi, embl, molecular, genetics, software, bioinformatics databases, genomics, sequencing, protein, computational biology, nucleotide, fasta, blast, srs, clustalw, bioinformatics research' /><meta name='author' content='EBI Web Services' /><meta name='Description' content='IMGT/HLA Download Page - FTP help page for downloading HLA sequences from the IMGT/HLA Database' /><meta http-equiv='Content-Language' content='en-GB' /><meta http-equiv='Window-target' content='_top' /><meta name='no-email-collection' content='http://www.unspam.com/noemailcollection/' /><meta name='generator' content='Dreamweaver 8' /><!-- InstanceBeginEditable name='doctitle' --><title>IMGT/HLA Database</title><!-- InstanceEndEditable --><link rel='stylesheet'  href='http://www.ebi.ac.uk/inc/css/contents.css'     type='text/css' /><link rel='stylesheet'  href='http://www.ebi.ac.uk/inc/css/userstyles.css'   type='text/css' /><script  src='http://www.ebi.ac.uk/inc/js/contents.js' type='text/javascript'></script><link rel='SHORTCUT ICON' href='http://www.ebi.ac.uk/bookmark.ico' /><!-- InstanceBeginEditable name='head' --><!--  start meta tags, css , javascript here   --> <style>code     {color: #6EA65A;         }div.code {text-align: left;          font-family: 'Courier New', Courier, mono ;          color: #616161;           padding: 12px;          background: #EEF6EC;          border: 1px solid #d9dadc;          margin: 0;}</style><script>function foo() {  alert('foo');}</script><style type='text/css'><!--.allelelink {font-family: 'Courier New', Courier, mono;	     font-size: 14px;	     color: #006666;}.normal     {font-family: 'Courier New', Courier, mono;	     font-size: 14px;	     color: #000000;}.spliced    {font-family: 'Courier New', Courier, mono;	     font-size: 14px;	     color: #000000;	     background: #edf5ea;}--></style> <!--  end meta tags, css , javascript here  --><!-- InstanceEndEditable --></head><body onload='if(navigator.userAgent.indexOf('MSIE') != -1) {document.getElementById('head').allowTransparency = true;}'>		<iframe src='/inc/head.html' name='head' id='head' frameborder='0' marginwidth='0px' marginheight='0px' scrolling='no'  width='100%' style='position:absolute; z-index: 1; height: 57px;'></iframe>	</div>	<div class='contents' id='contents'>			<table class='contentspane' id='contentspane' summary='The main content pane of the page' style='width: 100%'>				<tr>				  <td class='leftmargin'><img src='http://www.ebi.ac.uk/inc/images/spacer.gif' class='spacer' alt='spacer'  /></td>				  <td class='leftmenucell' id='leftmenucell'>				  	<div class='leftmenu' id='leftmenu'><img src='http://www.ebi.ac.uk/inc/images/spacer.gif' class='spacer' alt='spacer'  /></div>				  </td>				  <td class='contentsarea' id='contentsarea'>					<!-- InstanceBeginEditable name='contents' -->                    <!-- start contents here -->		<!-- start contents here -->				<div class='breadcrumbs'>                <a href='http://www.ebi.ac.uk/' class='firstbreadcrumb'>EBI</a>                <a href='/Databases/'>Databases</a>                <a href='/Databases/nucleotide.html'>Nucleotide Databases</a>                <a href='/imgt/hla/'>IMGT/HLA</a>                 <a href='/imgt/hla/align.html'>Sequence Alignments</a>               </div><h1>IMGT/HLA Database</h1>        <h2>Sequence Alignments based on Release 2.15.0 (06-October-2006)</h1>        <table border=0 width=100%>        <tr>        <td nowrap><ul></ul><span class='normal'><BR>&nbsp;<BR>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;10 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;20 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;30 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;40 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;50 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;60 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;70 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;80 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;90 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;100 <BR>&nbsp;A*01010101&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;GSHSMRYFFT SVSRPGRGEP RFIAVGYVDD TQFVRFDSDA ASQKMEPRAP WIEQEGPEYW DQETRNMKAH SQTDRANLGT LRGYYNQSED GSHTIQIMYG <BR>&nbsp;A*02010101&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;---------- ---------- ---------- ---------- ---R------ ---------- -G---KV--- ---H-VD--- ---------A ----V-R--- <BR>&nbsp;A*02010102L&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;---------- ---------- ---------- ---------- ---R------ ---------- -G---KV--- ---H-VD--- ---------A ----V-R--- <BR>&nbsp;A*020102&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;*--------- ---------- ---------- ---------- ---R------ ---------- -G---KV--- ---H-VD--- ---------A ----V-R--- <BR>&nbsp;A*020103&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;---------- ---------- ---------- ---------- ---R------ ---------- -G---KV--- ---H-VD--- ---------A ----V-R--- <BR>&nbsp;A*020104&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;*--------- ---------- ---------- ---------- ---R------ ---------- -G---KV--- ---H-VD--- ---------A ----V-R--- <BR>&nbsp;A*020105&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;*--------- ---------- ---------- ---------- ---R------ ---------- -G---KV--- ---H-VD--- ---------A ----V-R--- <BR>&nbsp;A*020106&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;*--------- ---------- ---------- ---------- ---R------ ---------- -G---KV--- ---H-VD--- ---------A ----V-R--- <BR>&nbsp;A*020107&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;*--------- ---------- ---------- ---------- ---R------ ---------- -G---KV--- ---H-VD--- ---------A ----V-R--- <BR>&nbsp;A*020108&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;*--------- ---------- ---------- ---------- ---R------ ---------- -G---KV--- ---H-VD--- ---------A ----V-R--- <BR>&nbsp;A*020109&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;*--------- ---------- ---------- ---------- ---R------ ---------- -G---KV--- ---H-VD--- ---------A ----V-R--- <BR>&nbsp;A*020110&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;*--------- ---------- ---------- ---------- ---R------ ---------- -G---KV--- ---H-VD--- ---------A ----V-R--- <BR>&nbsp;A*020111&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;*--------- ---------- ---------- ---------- ---R------ ---------- -G---KV--- ---H-VD--- ---------A ----V-R--- <BR>&nbsp;A*020112&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;*--------- ---------- ---------- ---------- ---R------ ---------- -G---KV--- ---H-VD--- ---------A ----V-R--- <BR>&nbsp;A*0205&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;--------Y- ---------- ---------- ---------- --RR------ ---------- -G---KV--- ---H-VD--- ---------A ----L-R--- <BR>&nbsp;A*020601&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;--------Y- ---------- ---------- ---------- ---R------ ---------- -G---KV--- ---H-VD--- ---------A ----V-R--- <BR>&nbsp;A*020602&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;*-------Y- ---------- ---------- ---------- ---R------ ---------- -G---KV--- ---H-VD--- ---------A ----V-R--- <BR>&nbsp;A*020603&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;*-------Y- ---------- ---------- ---------- ---R------ ---------- -G---KV--- ---H-VD--- ---------A ----V-R--- <BR>&nbsp;A*0215N&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;---------- ---------- ---------- ---------- ---R------ ---------- -G---KV--- ---H-VD--- ---------A ----V-R-C- <BR>&nbsp;A*03010101&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;---------- ---------- ---------- ---------- ---R------ ---------- ------V--Q -----VD--- ---------A ---------- <BR>&nbsp;A*03010102N&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;---------- ---------- ---------- ---------- ---R------ ---------- ------V--Q -----VD--- ---------A ---------- <BR>&nbsp;A*03010103&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;---------- ---------- ---------- ---------- ---R------ ---------- ------V--Q -----VD--- ---------A ---------- <BR>&nbsp;A*030102&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;*--------- ---------- ---------- ---------- ---R------ ---------- ------V--Q -----VD--- ---------A ---------- <BR>&nbsp;A*030103&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;*--------- ---------- ---------- ---------- ---R------ ---------- ------V--Q -----VD--- ---------A ---------- <BR>&nbsp;A*030104&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;*--------- ---------- ---------- ---------- ---R------ ---------- ------V--Q -----VD--- ---------A ---------- <BR>&nbsp;A*030105&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;*--------- ---------- ---------- ---------- ---R------ ---------- ------V--Q -----VD--- ---------A ---------- <BR>&nbsp;A*0303N&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;---------- ---------- ---------- ---------- ---R------ ---------- ------V--Q -----VD--- ---------A ---------- <BR>&nbsp;A*24020101&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;--------S- ---------- ---------- ---------- ---R------ ---------- -E--GKV--- -----E--RI ALR------A ----L-M-F- <BR>&nbsp;A*24020102L&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;--------S- ---------- ---------- ---------- ---R------ ---------- -E--GKV--- -----E--RI ALR------A ----L-M-F- <BR>&nbsp;A*240202&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;*-------S- ---------- ---------- ---------- ---R------ ---------- -E--GKV--- -----E--RI ALR------A ----L-M-F- <BR>&nbsp;A*240203&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;--------S- ---------- ---------- ---------- ---R------ ---------- -E--GKV--- -----E--RI ALR------A ----L-M-F- <BR>&nbsp;A*240204&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;--------S- ---------- ---------- ---------- ---R------ ---------- -E--GKV--- -----E--RI ALR------A ----L-M-F- <BR>&nbsp;A*240205&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;*-------S- ---------- ---------- ---------- ---R------ ---------- -E--GKV--- -----E--RI ALR------A ----L-M-F- <BR>&nbsp;A*240206&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;*-------S- ---------- ---------- ---------- ---R------ ---------- -E--GKV--- -----E--RI ALR------A ----L-M-F- <BR>&nbsp;A*240207&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;*-------S- ---------- ---------- ---------- ---R------ ---------- -E--GKV--- -----E--RI ALR------A ----L-M-F- <BR>&nbsp;A*240208&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;*-------S- ---------- ---------- ---------- ---R------ ---------- -E--GKV--- -----E--RI ALR------A ----L-M-F- <BR>&nbsp;A*240209&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;*-------S- ---------- ---------- ---------- ---R------ ---------- -E--GKV--- -----E--RI ALR------A ----L-M-F- <BR>&nbsp;A*240210&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;--------S- ---------- ---------- ---------- ---R------ ---------- -E--GKV--- -----E--RI ALR------A ----L-M-F- <BR>&nbsp;A*240211&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;*-------S- ---------- ---------- ---------- ---R------ ---------- -E--GKV--- -----E--RI ALR------A ----L-M-F- <BR>&nbsp;A*240212&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;*-------S- ---------- ---------- ---------- ---R------ ---------- -E--GKV--- -----E--RI ALR------A ----L-M-F- <BR>&nbsp;A*29010101&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;--------T- ---------- ---------- ---------- ---R------ ---------- -LQ---V--Q ---------- ---------A ------M--- <BR>&nbsp;A*29010102N&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;--------T- ---------- ---------- ---------- ---R------ ---------- -LQ---V--Q ---------- ---------A ------M--- <BR>&nbsp;A*300101&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;--------S- ------S--- ---------- ---------- ---R------ -----R---- ------V--Q -----VD--- ---------A ---------- <BR>&nbsp;A*300102&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;--------S- ------S--- ---------- ---------- ---R------ -----R---- ------V--Q -----VD--- ---------A ---------- <BR>&nbsp;A*4301&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;--------Y- ---------- ---------- ---------- ---R------ ---------- -LQ---V--- ---------- ---------- ------R--- <BR>&nbsp;A*680101&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;--------Y- ---------- ---------- ---------- ---R------ ---------- -RN---V--Q -----VD--- ---------A ------M--- <BR>&nbsp;A*680102&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;--------Y- ---------- ---------- ---------- ---R------ ---------- -RN---V--Q -----VD--- ---------A ------M--- <BR>&nbsp;A*680103&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;*-------Y- ---------- ---------- ---------- ---R------ ---------- -RN---V--Q -----VD--- ---------A ------M--- <BR>&nbsp;A*680104&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;*-------Y- ---------- ---------- ---------- ---R------ ---------- -RN---V--Q -----VD--- ---------A ------M--- <BR>&nbsp;A*7412N&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;*--------- ---------- ---------- ---------- ---R------ ---------- ------V--- -----VD--- ---------A ----..M--- <BR>&nbsp;A*8001&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;---------- ---------- ---------- S---Q----- ---R------ -----E---- -E----V--- ---N------ ---------- ---------- <BR>&nbsp;<BR>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;110 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;120 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;130 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;140 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;150 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;160 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;170 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;180 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;190 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;200 <BR>&nbsp;A*01010101&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;CDVGPDGRFL RGYRQDAYDG KDYIALNEDL RSWTAADMAA QITKRKWEAV HAAEQRRVYL EGRCVDGLRR YLENGKETLQ RTDPPKTHMT HHPISDHEAT <BR>&nbsp;A*02010101&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;----S-W--- ---H-Y---- ------K--- ---------- -T--H----A -V---L-A-- --T--EW--- ---------- ---A------ --AV------ <BR>&nbsp;A*02010102L&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;----S-W--- ---H-Y---- ------K--- ---------- -T--H----A -V---L-A-- --T--EW--- ---------- ---A------ --AV------ <BR>&nbsp;A*020102&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;----S-W--- ---H-Y---- ------K--- ---------- -T--H----A -V---L-A-- --T--EW--- ---------- --******** ********** <BR>&nbsp;A*020103&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;----S-W--- ---H-Y---- ------K--- ---------- -T--H----A -V---L-A-- --T--EW--- ---------- ---A------ --AV------ <BR>&nbsp;A*020104&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;----S-W--- ---H-Y---- ------K--- ---------- -T--H----A -V---L-A-- --T--EW--- ---------- --******** ********** <BR>&nbsp;A*020105&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;----S-W--- ---H-Y---- ------K--- ---------- -T--H----A -V---L-A-- --T--EW--- ---------- ---A------ --AV------ <BR>&nbsp;A*020106&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;----S-W--- ---H-Y---- ------K--- ---------- -T--H----A -V---L-A-- --T--EW--- ---------- ---A------ --AV------ <BR>&nbsp;A*020107&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;----S-W--- ---H-Y---- ------K--- ---------- -T--H----A -V---L-A-- --T--EW--- ---------- ---A------ --AV------ <BR>&nbsp;A*020108&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;----S-W--- ---H-Y---- ------K--- ---------- -T--H----A -V---L-A-- --T--EW--- ---------- ---A------ --AV------ <BR>&nbsp;A*020109&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;----S-W--- ---H-Y---- ------K--- ---------- -T--H----A -V---L-A-- --T--EW--- ---------- ---A------ --AV------ <BR>&nbsp;A*020110&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;----S-W--- ---H-Y---- ------K--- ---------- -T--H----A -V---L-A-- --T--EW--- ---------- ---A------ --AV------ <BR>&nbsp;A*020111&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;----S-W--- ---H-Y---- ------K--- ---------- -T--H----A -V---L-A-- --T--EW--- ---------- ---A------ --AV------ <BR>&nbsp;A*020112&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;----S-W--- ---H-Y---- ------K--- ---------- -T--H----A -V---L-A-- --T--EW--- ---------- --******** ********** <BR>&nbsp;A*0205&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;----S-W--- ---H-Y---- ------K--- ---------- -T--H----A -V---W-A-- --T--EW--- ---------- ---A------ --AV------ <BR>&nbsp;A*020601&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;----S-W--- ---H-Y---- ------K--- ---------- -T--H----A -V---L-A-- --T--EW--- ---------- ---A------ --AV------ <BR>&nbsp;A*020602&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;----S-W--- ---H-Y---- ------K--- ---------- -T--H----A -V---L-A-- --T--EW--- ---------- ---A------ --AV------ <BR>&nbsp;A*020603&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;----S-W--- ---H-Y---- ------K--- ---------- -T--H----A -V---L-A-- --T--EW--- ---------- ---A------ --AV------ <BR>&nbsp;A*0215N&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;----S-W--- ---H-Y---- ------K--- ---------- -T--H----A -V---L-A-- --T--EW--- ---------- ---A------ --AV------ <BR>&nbsp;A*03010101&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;----S----- ---------- ---------- ---------- ---------A -E---L-A-- D-T--EW--- ---------- ---------- ---------- <BR>&nbsp;A*03010102N&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;----S----- ---------- ---------- ---------- ---------A -E---L-A-- D-T--EW--- ---------- ---------- ---------- <BR>&nbsp;A*03010103&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;----S----- ---------- ---------- ---------- ---------A -E---L-A-- D-T--EW--- ---------- ---------- ---------- <BR>&nbsp;A*030102&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;----S----- ---------- ---------- ---------- ---------A -E---L-A-- D-T--EW--- ---------- --******** ********** <BR>&nbsp;A*030103&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;----S----- ---------- ---------- ---------- ---------A -E---L-A-- D-T--EW--- ---------- --******** ********** <BR>&nbsp;A*030104&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;----S----- ---------- ---------- ---------- ---------A -E---L-A-- D-T--EW--- ---------- --******** ********** <BR>&nbsp;A*030105&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;----S----- ---------- ---------- ---------- ---------A -E---L-A-- D-T--EW--- ---------- --******** ********** <BR>&nbsp;A*0303N&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;..--S----- ---------- ---------- ---------- ---------A -E---L-A-- D-T--EW--- ---------- ---------- ---------- <BR>&nbsp;A*24020101&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;----S----- ---H-Y---- ------K--- ---------- ---------A -V---Q-A-- --T------- ---------- ---------- ---------- <BR>&nbsp;A*24020102L&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;----S----- ---H-Y---- ------K--- ---------- ---------A -V---Q-A-- --T------- ---------- ---------- ---------- <BR>&nbsp;A*240202&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;----S----- ---H-Y---- ------K--- ---------- ---------A -V---Q-A-- --T------- ---------- --******** ********** <BR>&nbsp;A*240203&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;----S----- ---H-Y---- ------K--- ---------- ---------A -V---Q-A-- --T------- ---------- ---------- ---------- <BR>&nbsp;A*240204&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;----S----- ---H-Y---- ------K--- ---------- ---------A -V---Q-A-- --T------- ---------- ---------- ---------- <BR>&nbsp;A*240205&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;----S----- ---H-Y---- ------K--- ---------- ---------A -V---Q-A-- --T------- ---------- --******** ********** <BR>&nbsp;A*240206&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;----S----- ---H-Y---- ------K--- ---------- ---------A -V---Q-A-- --T------- ---------- ---------- ---------- <BR>&nbsp;A*240207&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;----S----- ---H-Y---- ------K--- ---------- ---------A -V---Q-A-- --T------- ---------- ---------- ---------- <BR>&nbsp;A*240208&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;----S----- ---H-Y---- ------K--- ---------- ---------A -V---Q-A-- --T------- ---------- --******** ********** <BR>&nbsp;A*240209&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;----S----- ---H-Y---- ------K--- ---------- ---------A -V---Q-A-- --T------- ---------- --******** ********** <BR>&nbsp;A*240210&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;----S----- ---H-Y---- ------K--- ---------- ---------A -V---Q-A-- --T------- ---------- ---------- ---------- <BR>&nbsp;A*240211&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;----S----- ---H-Y---- ------K--- ---------- ---------A -V---Q-A-- --T------- ---------- --******** ********** <BR>&nbsp;A*240212&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;----S----- ---H-Y---- ------K--- ---------- ---------A -V---Q-A-- --T------- ---------- --******** ********** <BR>&nbsp;A*29010101&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;-H--S----- ---------- ---------- ---------- ---Q-----A RV---L-A-- --T--EW--- ---------- ---A------ --AV------ <BR>&nbsp;A*29010102N&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;-H--S----- ---------- ---------- ---------- ---Q-----A RV---L-A-- --T--EW--- ---------- ---A------ --AV------ <BR>&nbsp;A*300101&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;----S----- ---E-H---- ---------- ---------- ---Q-----A RW---L-A-- --T--EW--- ---------- ---------- ---------- <BR>&nbsp;A*300102&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;----S----- ---E-H---- ---------- ---------- ---Q-----A RW---L-A-- --T--EW--- ---------- ---------- ---------- <BR>&nbsp;A*4301&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;---------- ---Q------ ---------- ---------- ---Q----TA -E---W-A-- -----EW--- ---------- ---A------ --AV------ <BR>&nbsp;A*680101&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;----S----- ---------- ------K--- ---------- -T--H----A -V---W-A-- --T--EW--- ---------- ---A------ --AV------ <BR>&nbsp;A*680102&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;----S----- ---------- ------K--- ---------- -T--H----A -V---W-A-- --T--EW--- ---------- ---A------ --AV------ <BR>&nbsp;A*680103&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;----S----- ---------- ------K--- ---------- -T--H----A -V---W-A-- --T--EW--- ---------- --******** ********** <BR>&nbsp;A*680104&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;----S----- ---------- ------K--- ---------- -T--H----A -V---W-A-- --T--EW--- ---------- --******** ********** <BR>&nbsp;A*7412N&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;--------L- ---Q------ ---------- ---------- ---Q-----A RV---L-A-- --T--EW--- ---------- --******** ********** <BR>&nbsp;A*8001&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;----S----- ---------- ---------- ---------- ---------A RR---L-A-- --E------- ---------- ---------- ---------- <BR>&nbsp;<BR>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;210 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;220 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;230 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;240 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;250 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;260 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;270 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;280 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;290 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;300 <BR>&nbsp;A*01010101&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;LRCWALGFYP AEITLTWQRD GEDQTQDTEL VETRPAGDGT FQKWAAVVVP SGEEQRYTCH VQHEGLPKPL TLRWELSSQP TIPIVGIIAG LVLLGAVITG <BR>&nbsp;A*02010101&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;------S--- ---------- ---------- ---------- ---------- --Q------- ---------- -----P---- ---------- ---F------ <BR>&nbsp;A*02010102L&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;------S--- ---------- ---------- ---------- ---------- --Q------- ---------- -----P---- ---------- ---F------ <BR>&nbsp;A*020102&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;********** ********** ********** ********** ********** ********** ********** ********** ********** ********** <BR>&nbsp;A*020103&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;------S--- ---------- ---------- ---------- ---------- --Q------- ---------- ----****** ********** ********** <BR>&nbsp;A*020104&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;********** ********** ********** ********** ********** ********** ********** ********** ********** ********** <BR>&nbsp;A*020105&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;------S--- ---------- ---------- ---------- ---------- --Q------- ---------- ----****** ********** ********** <BR>&nbsp;A*020106&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;------S--- ---------- ---------- ---------- ---------- --Q------- ---------- ----****** ********** ********** <BR>&nbsp;A*020107&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;------S--- ---------- ---------- ---------- ---------- --Q------- ---------- ----****** ********** ********** <BR>&nbsp;A*020108&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;------S--- ---------- ---------- ---------- ---------- --Q------- ---------- ----****** ********** ********** <BR>&nbsp;A*020109&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;------S--- ---------- ---------- ---------- ---------- --Q------- ---------- ----****** ********** ********** <BR>&nbsp;A*020110&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;------S--- ---------- ---------- ---------- ---------- --Q------- ---------- ----****** ********** ********** <BR>&nbsp;A*020111&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;------S--- ---------- ---------- ---------- ---------- --Q------- ---------- ----****** ********** ********** <BR>&nbsp;A*020112&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;********** ********** ********** ********** ********** ********** ********** ********** ********** ********** <BR>&nbsp;A*0205&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;------S--- ---------- ---------- ---------- ---------- --Q------- ---------- -----P---- ---------- ---F------ <BR>&nbsp;A*020601&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;------S--- ---------- ---------- ---------- ---------- --Q------- ---------- -----P---- ---------- ---F------ <BR>&nbsp;A*020602&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;------S--- ---------- ---------- ---------- ---------- --Q------- ---------- ----****** ********** ********** <BR>&nbsp;A*020603&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;------S--- ---------- ---------- ---------- ---------- --Q------- ---------- ----****** ********** ********** <BR>&nbsp;A*0215N&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;------S--- ---------- ---------- ---------- ---------- --Q---X     <BR>&nbsp;A*03010101&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;---------- ---------- ---------- ---------- ---------- ---------- ---------- ---------- ---------- ---------- <BR>&nbsp;A*03010102N&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;---------- ---------- ---------- ---------- ---------- ---------- ---------- ----DKEGDG GVMSLRESRS -SGDLX <BR>&nbsp;A*03010103&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;---------- ---------- ---------- ---------- ---------- ---------- ---------- ---------- ---------- ---------- <BR>&nbsp;A*030102&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;********** ********** ********** ********** ********** ********** ********** ********** ********** ********** <BR>&nbsp;A*030103&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;********** ********** ********** ********** ********** ********** ********** ********** ********** ********** <BR>&nbsp;A*030104&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;********** ********** ********** ********** ********** ********** ********** ********** ********** ********** <BR>&nbsp;A*030105&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;********** ********** ********** ********** ********** ********** ********** ********** ********** ********** <BR>&nbsp;A*0303N&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;---------- ---------- ---------- ---------- ---------- ---------- ---------- ---------- ---------- ---------- <BR>&nbsp;A*24020101&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;---------- ---------- ---------- ---------- ---------- ---------- ---------- -----P---- -V-------- ---------- <BR>&nbsp;A*24020102L&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;---------- ---------- ---------- ---------- ---------- ---------- ---------- -----P---- -V-------- ---------- <BR>&nbsp;A*240202&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;********** ********** ********** ********** ********** ********** ********** ********** ********** ********** <BR>&nbsp;A*240203&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;---------- ---------- ---------- ---------- ---------- ---------- ---------- -----P---- -V-------- ---------- <BR>&nbsp;A*240204&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;---------- ---------- ---------- ---------- ---------- ---------- ---------- ----****** ********** ********** <BR>&nbsp;A*240205&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;********** ********** ********** ********** ********** ********** ********** ********** ********** ********** <BR>&nbsp;A*240206&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;---------- ---------- ---------- ---------- ---------- ---------- ---------- ----****** ********** ********** <BR>&nbsp;A*240207&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;---------- ---------- ---------- ---------- ---------- ---------- ---------- ----****** ********** ********** <BR>&nbsp;A*240208&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;********** ********** ********** ********** ********** ********** ********** ********** ********** ********** <BR>&nbsp;A*240209&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;********** ********** ********** ********** ********** ********** ********** ********** ********** ********** <BR>&nbsp;A*240210&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;---------- ---------- ---------- ---------- ---------- ---------- ---------- -----P---- -V-------- ---------- <BR>&nbsp;A*240211&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;********** ********** ********** ********** ********** ********** ********** ********** ********** ********** <BR>&nbsp;A*240212&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;********** ********** ********** ********** ********** ********** ********** ********** ********** ********** <BR>&nbsp;A*29010101&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;------S--- ---------- ---------- ---------- -----S---- --Q------- ---------- -----P---- ---------- ---F---FA- <BR>&nbsp;A*29010102N&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;------S--- ---------- ---------- ---------- -----S---- --Q------- ---------- ----VKEGDG GVMSFRESRS -S.DLX <BR>&nbsp;A*300101&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;---------- ---------- ---------- ---------- ---------- ---------- ---------- ---------- ---------- ---------- <BR>&nbsp;A*300102&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;---------- ---------- ---------- ---------- ---------- ---------- ---------- ----****** ********** ********** <BR>&nbsp;A*4301&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;------S--- ---------- ---------- ---------- -----S---- --Q------- ---------- -----P---- ---------- ---F----A- <BR>&nbsp;A*680101&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;------S--- ---------- ---------- ---------- ----V----- --Q------- ---------- -----P---- ---------- ---F------ <BR>&nbsp;A*680102&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;------S--- ---------- ---------- ---------- ----V----- --Q------- ---------- -----P---- ---------- ---F------ <BR>&nbsp;A*680103&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;********** ********** ********** ********** ********** ********** ********** ********** ********** ********** <BR>&nbsp;A*680104&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;********** ********** ********** ********** ********** ********** ********** ********** ********** ********** <BR>&nbsp;A*7412N&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;********** ********** ********** ********** ********** ********** ********** ********** ********** ********** <BR>&nbsp;A*8001&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;------S--- ---------- ---------- ---------- ---------- --K-K----- -------E-- -----P---- ---------- --------A- <BR>&nbsp;<BR>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;310 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;320 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;330 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;340 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;350 <BR>&nbsp;A*01010101&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;AVVAAVMWRR KSSDRKGGSY TQAASSDSAQ GSDVSLTACK V&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;      <BR>&nbsp;A*02010101&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;---------- ---------- S--------- ---------- -&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;      <BR>&nbsp;A*02010102L&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;---------- ---------- S--------- ---------- -&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;      <BR>&nbsp;A*020102&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;********** ********** ********** ********** *&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;      <BR>&nbsp;A*020103&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;********** ********** ********** ********** *&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;      <BR>&nbsp;A*020104&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;********** ********** ********** ********** *&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;      <BR>&nbsp;A*020105&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;********** ********** ********** ********** *&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;      <BR>&nbsp;A*020106&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;********** ********** ********** ********** *&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;      <BR>&nbsp;A*020107&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;********** ********** ********** ********** *&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;      <BR>&nbsp;A*020108&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;********** ********** ********** ********** *&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;      <BR>&nbsp;A*020109&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;********** ********** ********** ********** *&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;      <BR>&nbsp;A*020110&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;********** ********** ********** ********** *&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;      <BR>&nbsp;A*020111&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;********** ********** ********** ********** *&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;      <BR>&nbsp;A*020112&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;********** ********** ********** ********** *&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;      <BR>&nbsp;A*0205&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;---------- ---------- S--------- ---------- -&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;      <BR>&nbsp;A*020601&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;---------- ---------- S--------- ---------- -&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;      <BR>&nbsp;A*020602&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;********** ********** ********** ********** *&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;      <BR>&nbsp;A*020603&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;********** ********** ********** ********** *&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;      <BR>&nbsp;A*0215N&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;    &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;      <BR>&nbsp;A*03010101&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;---------- ---------- ---------- ---------- -&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;      <BR>&nbsp;A*03010102N&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;    &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;      <BR>&nbsp;A*03010103&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;---------- ---------- ---------- ---------- -&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;      <BR>&nbsp;A*030102&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;********** ********** ********** ********** *&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;      <BR>&nbsp;A*030103&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;********** ********** ********** ********** *&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;      <BR>&nbsp;A*030104&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;********** ********** ********** ********** *&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;      <BR>&nbsp;A*030105&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;********** ********** ********** ********** *&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;      <BR>&nbsp;A*0303N&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;---------- ---------- ---------- ---------- -&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;      <BR>&nbsp;A*24020101&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;---------- N--------- S--------- ---------- -&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;      <BR>&nbsp;A*24020102L&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;---------- N--------- S--------- ---------- -&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;      <BR>&nbsp;A*240202&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;********** ********** ********** ********** *&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;      <BR>&nbsp;A*240203&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;---------- N--------- S--------- ---------- -&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;      <BR>&nbsp;A*240204&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;********** ********** ********** ********** *&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;      <BR>&nbsp;A*240205&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;********** ********** ********** ********** *&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;      <BR>&nbsp;A*240206&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;********** ********** ********** ********** *&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;      <BR>&nbsp;A*240207&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;********** ********** ********** ********** *&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;      <BR>&nbsp;A*240208&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;********** ********** ********** ********** *&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;      <BR>&nbsp;A*240209&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;********** ********** ********** ********** *&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;      <BR>&nbsp;A*240210&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;---------- N--------- S--------- ---------- -&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;      <BR>&nbsp;A*240211&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;********** ********** ********** ********** *&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;      <BR>&nbsp;A*240212&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;********** ********** ********** ********** *&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;      <BR>&nbsp;A*29010101&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;------R--- ---------- S--------- ---M------ -&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;      <BR>&nbsp;A*29010102N&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;    &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;      <BR>&nbsp;A*300101&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;---------- ---------- ---------- ---------- -&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;      <BR>&nbsp;A*300102&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;********** ********** ********** ********** *&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;      <BR>&nbsp;A*4301&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;---------- ---------- S--------- ---M------ -&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;      <BR>&nbsp;A*680101&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;---------- ---------- S--------- ---------- -&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;      <BR>&nbsp;A*680102&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;---------- ---------- S--------- ---------- -&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;      <BR>&nbsp;A*680103&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;********** ********** ********** ********** *&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;      <BR>&nbsp;A*680104&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;********** ********** ********** ********** *&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;      <BR>&nbsp;A*7412N&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;********** ********** ********** ********** *&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;      <BR>&nbsp;A*8001&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;---------K ---V------ S--------- ---------- -&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;      </span></td></tr></table>       <h2>Further Information</h2>       <p class='textjustify'>For more information about the database or for all queries please contact       <a href='/support/ipd.php' class='greenbold'>IPD Support</a>.</p><!-- InstanceEndEditable -->										<img src='http://www.ebi.ac.uk/inc/images/spacer.gif' class='spacer'  alt='spacer' /></td>					<td class='rightmenucell' id='rightmenucell'>					  <div class='rightmenu' id='rightmenu'><img src='http://www.ebi.ac.uk/inc/images/spacer.gif' class='spacer' alt='spacer' /></div>				  	</td>				</tr>				</table>				<table class='footerpane' id='footerpane' summary='The main footer pane of the page'>				<tr>				  <td colspan ='4' class='footerrow'>					<div class='footerdiv' id='footerdiv'  style='z-index:2;'>						<iframe src='/inc/foot.html' name='foot' frameborder='0' marginwidth='0px' marginheight='0px' scrolling='no'  height='22px' width='100%'  style='z-index:2;'></iframe>					</div>				  </td>				</tr>	  </table>	  </table>	  <script  src='http://www.ebi.ac.uk/inc/js/footer.js' type='text/javascript'></script>	</div></body><!-- InstanceEnd --></html>";

        //ONLINE doesnt work.
		if (ONLINE){
			//Create connection
			URL url;
			final String inUrl = "http://www.ebi.ac.uk/cgi-bin/imgt/hla/align.cgi";
			try {
				url = new URL(inUrl);
				URLConnection connection = url.openConnection();
				HttpURLConnection httpConn = (HttpURLConnection) connection;
				
				//Convert allele array to a POST message
				StringBuffer alleleStr = new StringBuffer();
				for (int i = 0; i<allele_count; i++){    
					alleleStr.append(allelesAsc.get(i));
					if (i<allele_count - 1)
						alleleStr.append("%2C");
				}
				
				if (!DEBUG){
					//Build POST string
                    String msg = "Type=MatureProtein&submit=Align+Sequences+Now&Sequences="+alleleStr.toString()+"&Reference="+ref+"&Printing=P&Omit=N&gene="+ref+"&Formatting=10&Display=Show+Mismatches";
                    byte[]b = msg.getBytes();
					
					//Set HTTP parameters.
					httpConn.setRequestProperty("Content-Length", String.valueOf(b.length));
					httpConn.setRequestProperty("Content-Type","application/x-www-form-urlencoded");
					httpConn.setRequestMethod("POST");
					httpConn.setDoOutput(true);
					httpConn.setDoInput(true);
					
					//Send message
					OutputStream out = httpConn.getOutputStream();
					out.write(b);    
					out.close();
					
					//Read response
					InputStreamReader isr = new InputStreamReader(httpConn.getInputStream());
					BufferedReader in = new BufferedReader(isr);
					StringBuffer sb = new StringBuffer();
					
					String inputLine;
					while ((inputLine = in.readLine()) != null)
						sb.append(inputLine);	
					in.close();
					response = sb.toString();
				}
			} catch (Exception e) {
				if (!Settings.SUPPRESS_WARNINGS)Tools.printErr("\nUnable to establish an internet connection. Please make sure your internet connection is active and try again.");
				return;
			}
		}


        // PJ - DOWNLOAD is the only thing that works..
        if(DOWNLOAD){

            FileDownloader fd = new FileDownloader();
            String locus = String.valueOf(loc);
            //try to see if the file exists "<locus>.txt". else remove the last character from locus and try again.. i.e if locus is DRB1, then it checks for DRB1 and then for DRB (eliminating the "1")
            try{
                if(EnumUtils.isValidEnum(Locus.class, String.valueOf(loc))){
                    fd.setExternalFile("ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/alignments/"+locus+"_prot.txt");
                    fd.setLocalFile(locus + "_prot.txt");
                    fd.setMaxRetryCount(1);
                    fd.download();
                }

            }catch (FileDownloader.PermanentDownloadErrorException e){
                try{
                    if(e.getStackTrace().length>0){
                        String newlocus = locus.substring(0, locus.length()-1);

                        fd.setExternalFile("ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/alignments/"+newlocus+"_prot.txt");
                        fd.setLocalFile(locus + "_prot.txt");
                        fd.setMaxRetryCount(1);
                        fd.download();

                    }
                }catch(FileDownloader.PermanentDownloadErrorException ex){
                    ex.printStackTrace();
                }
            }

            System.out.println("done with downloading");

            try{
                int indexOf1=0;
                String filename = locus + "_prot.txt";
                BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(filename)));
                StringBuffer sb = new StringBuffer();
                String inputLine;
                while ((inputLine = in.readLine()) != null){
                    // split each line by locus
                    if(inputLine.startsWith(" "+loc)){

                        //find index of position where the real protein starts at position 1.
                        if(inputLine.length()>indexOf1){
                            String pattern = loc+"\\*";
                            String inputLineaddingSplit = inputLine.substring(0, indexOf1-1)+"$"+inputLine.substring(indexOf1, inputLine.length());

                            String modifiedInputLine = inputLineaddingSplit.replaceAll(pattern,"");
                            modifiedInputLine = modifiedInputLine.substring(1, modifiedInputLine.length()-1);

                            sb.append(modifiedInputLine);
                            sb.append("!");
                        }else{
                            System.out.println("allele sequence length"+ inputLine.length()+" shorter than indexof1 position:"+indexOf1);
                        }
                    }else if(inputLine.startsWith(" Prot")){

                        String regExToFind = "[1-9]\\d*$";
                        //pj test
                        //indexOf1 = inputLine.indexOf("1");
                        int oldIndexOf1 = inputLine.indexOf("1");
                        indexOf1 = findIndexOf(Pattern.compile(regExToFind), inputLine);
                        System.out.println("the indexOf position is:"+indexOf1 +  "--versus the old indexof position is:" + oldIndexOf1);
                    }
                }
                in.close();
                modResponse = sb.toString()+"hello";
            }catch(Exception e){
                // file does not exists
                e.printStackTrace();
            }
        }

        //not used anymore..
		else{
			try{
				String filename = "align" + loc + ".htm";
	        	URL url = HLAAlign.class.getResource(filename);
	        	BufferedReader in = new BufferedReader(new InputStreamReader(url.openStream()));
	        	StringBuffer sb = new StringBuffer();
	        	String inputLine;
		    	while ((inputLine = in.readLine()) != null)
		    	    sb.append(inputLine);
		    	in.close();
	        	response = sb.toString();
	        }catch(Exception e){
	        	// file does not exists
	        	e.printStackTrace();
	        }
		}


        //create a linkedHashMap to store the retrieved alleles, the retrieved headers and the data from the locus file.


        Map<String, String> alleleMap = new LinkedHashMap<>();
        Map<String, String> headerMap = new LinkedHashMap<>();
        Map<String, String> locusFileAlleleMap = new LinkedHashMap<>();

        try{
            String[] alleleStr = modResponse.split("!", -1);

            for(String eachSeq : alleleStr){
                if(!eachSeq.equalsIgnoreCase("hello")){
                    String headerAndSeq = eachSeq.replaceFirst(" ", "&");
                    String header = headerAndSeq.split("&", -1)[0];
                    String[] headerFieldArray = header.split(":", -1);
                    String seq = headerAndSeq.split("\\&", -1)[1].split("\\$", -1)[1].replaceAll(" ","");


                    if(locusFileAlleleMap.size()>0){
                        if(locusFileAlleleMap.containsKey(header)){
                            String sequence = locusFileAlleleMap.get(header);
                            String s = sequence+seq;
                            locusFileAlleleMap.put(header, s);
                        }else{
                            locusFileAlleleMap.put(header, seq);
                        }
                    }else{
                        locusFileAlleleMap.put(header, seq);
                    }
                }
            }

            System.out.println("size of the locusFileMap hash is:"+locusFileAlleleMap.size());

            //set ref seq in headermap. remove allele objects that werent found in imgt locust.txt file.
            if(locusFileAlleleMap.containsKey(ref)){
                headerMap.put(ref, locusFileAlleleMap.get(ref));
            }

            //set alleleMap
            int alSerial=0;
            allelesTobeRemoved = new ArrayList<String>(allele_count);

            for(String allele : allelesAsc){

                int foundAllele=0;

                if(locusFileAlleleMap.containsKey(allele)){
                    alleleMap.put(allele, locusFileAlleleMap.get(allele));
                    foundAllele=1;
                }else{


                    //if 100% string match doesnt find allele.. chop allele down to two fields separated by : and check if two fields match... if not bring it down to one field and check if one field matches..
                    Set<String> collectionLocusheaders = locusFileAlleleMap.keySet();
                    Iterator<String> iter = collectionLocusheaders.iterator();

                    String[] alleleFieldArr = allele.split(":", -1);

                    //
                    if(alleleFieldArr.length>=2){

                        //if 100% string match isnt possible, then disregard the alleles ending in N when you clip the fields to 2 places.

                        while(iter.hasNext()){
                            String locusAll = iter.next();
                            if(!locusAll.endsWith("N")){
                                String[] locusAllFieldArr = locusAll.split(":", -1);
                                if(alleleFieldArr[0].equalsIgnoreCase(locusAllFieldArr[0]) &&
                                        alleleFieldArr[1].equalsIgnoreCase(locusAllFieldArr[1])){
                                    alleleMap.put(allele, locusFileAlleleMap.get(locusAll));
                                    foundAllele=1;
                                    break;
                                }
                            }
                        }

                    }else if(alleleFieldArr.length==1){
                        while(iter.hasNext()){
                            String locusAll = iter.next();
                            if(!locusAll.endsWith("N")){
                                String[] locusAllFieldArr = locusAll.split(":", -1);
                                if(alleleFieldArr[0].equalsIgnoreCase(locusAllFieldArr[0])){
                                    alleleMap.put(allele, locusFileAlleleMap.get(locusAll));
                                    foundAllele=1;
                                    break;
                                }
                            }

                        }
                    }


                }

                //remove alleles if allele not found.
                if(foundAllele==0){
                    System.out.println("could not find allele:"+allele);
                    if (!Settings.SUPPRESS_WARNINGS)Tools.printErr("\n*Allele '" + allelesAsc.indexOf(allele) + "' does not exist for HLA-" + loc + " and will be ignored");
                    deltas.remove(allelesIn.indexOf(allele));
                    allelesIn.remove(allelesIn.indexOf(allele));
                    allelesTobeRemoved.add(alSerial, allele);
                    allele_count=allele_count-1;
                }

                alSerial+=1;
            }


            for(int r=0; r<allelesTobeRemoved.size(); r++){
                allelesAsc.remove(allelesAsc.indexOf(allelesTobeRemoved.get(r)));
            }

            int serial=0;

            //commented out for now..
            Set<String> headerKeyset = headerMap.keySet();
            Iterator<String> iter = headerKeyset.iterator();
            while (iter.hasNext()) {
                String header = iter.next();
                //aligns[serial]=headerMap.get(header).toCharArray();
                //allelesRetrieved.add(serial, header);
                //serial+=1;
                refseq = headerMap.get(header).toCharArray();
                System.out.println("done with getting sequences, now alignment:"+ header);
            }


            //populate allelesRetrieved
            Set<String> alleleKeyset = alleleMap.keySet();
            Iterator<String> aliter = alleleKeyset.iterator();
            while (aliter.hasNext()) {
                String alheader = aliter.next();
                aligns[serial]=alleleMap.get(alheader).toCharArray();
                allelesRetrieved.add(serial, alheader);
                serial+=1;
                System.out.println("done with getting sequences, now alignment:"+ alheader);
            }

        }



		//Parse HLA alignment
        catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	private void printError(){
		if (!Settings.SUPPRESS_WARNINGS)Tools.printErr("\n\n***Locus HLA-" + loc + " has no alleles eligible for analysis!");
	}
	
	private String cleanHTML(String in){
		String ret;
		if (in == null || in.indexOf(loc+"*"+ref) == -1)
			return null;
        ret = in.substring(in.indexOf("&nbsp;"+loc+"*"));
        ret = ret.substring(0, ret.indexOf("</"));
        ret = ret.replaceAll("<BR>","!");
        ret = ret.replaceAll("<br>","!");
        ret = ret.replaceAll("&nbsp;"+loc+"\\*", "");
        ret = ret.replaceAll("&nbsp;", "");
        ret = ret.replaceAll("\\s", "");
        return ret;
	}

    //modified sortAligns
    private void sortAligns(){
        char[][] sortedAligns = new char[allele_count][aligns[0].length];
        if (deltas!=null){
            int deltaIndex;
            for (int i = 0; i < allelesRetrieved.size(); i++){
                deltaIndex = allelesAsc.indexOf(allelesIn.get(i));
                sortedAligns[i] = aligns[deltaIndex];
            }
            aligns = sortedAligns;
        }
    }

    /*private void sortAligns(){
		char[][] sortedAligns = new char[allele_count][aligns[0].length];
		if (deltas!=null){
			int deltaIndex;
			for (int i = 0; i < allelesRetrieved.size(); i++){
				deltaIndex = allelesAsc.indexOf(allelesIn.get(i));
				sortedAligns[i] = aligns[deltaIndex];
			}
			aligns = sortedAligns;
		}
	}*/
	
	private int alleleComp(String aligned, String query){
		if (aligned.equals(query))
			return 0;
		else{
			if (query.length() == 2){ //lores
				return -1;
			}
			else if (aligned.length() >= query.length()){
				if (aligned.startsWith(query))
					return 0;
				else if(aligned.substring(0, query.length()).compareTo(query) > 0)
					return 1;
				else
					return -1;
			}
			else if (aligned.length() < query.length()){
				return -1;
			}
		}
		return 0;
	}
	
	public String blanks(int n){
		StringBuffer ret = new StringBuffer();
		for (int i =0; i<n; i++)
			ret.append(" ");
		return ret.toString();
	}
	
	public String toString(){
		if (allelesRetrieved == null || allelesRetrieved.isEmpty()) return Settings.EMPTY_STRING;
		StringBuffer sb = new StringBuffer();
		String header;

        //this value changes based on length of ref seq
        ALIGNMENTS_PRINT_UNTIL_POS=refseq.length;
        System.out.println("alignment until pos is: "+ ALIGNMENTS_PRINT_UNTIL_POS);

        //changed : to = in the statement below.
		sb.append("\n\nAlignment of HLA-" + loc + " alleles=");
		sb.append("\n" + (deltas!=null?"Delta\t":"") + "Allele" + blanks(MAX_HEADER-6));
		for (int i = MAX_HEADER; i<=Math.min(ALIGNMENTS_PRINT_UNTIL_POS, COLS*COL); i=i+10){
			sb.append("    ." + blanks(5 - String.valueOf(i).length()-1) + i + "|");
		}
		if (deltas!=null){
			boolean flag = false;
			for (int i = 0; i < allelesRetrieved.size(); i++){
                header = allelesIn.get(i);
                //header = allelesRetrieved.get(i);
				if (!flag && deltas.get(i).charAt(0) == '-') {
                    sb.append("\n");
                    flag = true;
                }

                /*System.out.println(deltas.get(i));
                System.out.println(header);
                System.out.println(blanks(MAX_HEADER-header.length()));
                System.out.println(String.valueOf(aligns[i])+" \nand size is:"+ aligns[i].length);
                System.out.println(Math.min(ALIGNMENTS_PRINT_UNTIL_POS,COLS*COL));
                System.out.println(String.valueOf(aligns[i]).substring(0, Math.min(ALIGNMENTS_PRINT_UNTIL_POS, COLS * COL)));*/

				sb.append("\n" + deltas.get(i) + "\t" + header + blanks(MAX_HEADER - header.length()) + String.valueOf(aligns[i]).substring(0, Math.min(ALIGNMENTS_PRINT_UNTIL_POS, COLS * COL)));
                //System.out.println(sb.toString());
			}
		}
		else{
			for (int i = 0; i<allelesAsc.size(); i++){
				header = allelesAsc.get(i);
				sb.append("\n" + header + blanks(MAX_HEADER-header.length()) + String.valueOf(aligns[i]).substring(0, Math.min(ALIGNMENTS_PRINT_UNTIL_POS,COLS*COL)));
			}
		}
		return sb.toString();
	}




    //find index of this pattern ^[1-9]\d*$
    /** @return index of pattern in s or -1, if not found */
    public static int findIndexOf(Pattern pattern, String s) {
        Matcher matcher = pattern.matcher(s);
        return matcher.find() ? matcher.start() : -1;
    }



}
