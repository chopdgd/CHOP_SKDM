package utils;

public class FisherExact {

	/* This script can be used to  test statistically whether
	 * there is any relation between two categorical variables (with two levels). 
	 * Adopted from Oyvind Langsrud's javascript (http://www.matforsk.no/OLA/FISHER.HTM)
	 * by Stathis Kanterakis
	 * 
	 * Copyright (C) Oyvind Langsrud : All right reserved.
	 */
	private static double sn11,sn1_,sn_1,sn,sprob;
	private static double sleft,sright,sless,slarg;
	private static double left,right,twotail;
	private static int n11old=-1;
	private static int n12old=-1;
	private static int n21old=-1;
	private static int n22old=-1;
	
	public static final boolean DEBUG = false;
	
	private static double lngamm (double z){
		/* Reference: "Lanczos, C. 'A precision approximation
		 * of the gamma function', J. SIAM Numer. Anal., B, 1, 86-96, 1964."
		 * Translation of  Alan Miller's FORTRAN-implementation
		 * See http://lib.stat.cmu.edu/apstat/245
		 */
		double x = 0.0;
		x += 0.1659470187408462e-06/(z+7);
		x += 0.9934937113930748e-05/(z+6);
		x -= 0.1385710331296526    /(z+5);
		x += 12.50734324009056     /(z+4);
		x -= 176.6150291498386     /(z+3);
		x += 771.3234287757674     /(z+2);
		x -= 1259.139216722289     /(z+1);
		x += 676.5203681218835     /(z);
		x += 0.9999999999995183;
		return (Math.log(x)-5.58106146679532777-z+(z-0.5)*Math.log(z+6.5));
	}

	private static double lnfact(double n){
		if(n<=1) return(0.0);
		return(lngamm(n+1));
	}

	private static double lnbico(double n, double k){
		return(lnfact(n)-lnfact(k)-lnfact(n-k));
	}

	private static double hyper_323(double n11,double n1_,double n_1,double n){
		return(Math.exp(lnbico(n1_,n11)+lnbico(n-n1_,n_1-n11)-lnbico(n,n_1)));
	}
	
	private static double hyper(double n11){
		return(hyper0(n11,0.0,0.0,0.0));
	}

	private static double hyper0(double n11i,double n1_i,double n_1i,double ni){
		if(!(n1_i!=0|n_1i!=0|ni!=0)){
			if(!(n11i % 10 == 0)){
				if(n11i==sn11+1){
					sprob *= ((sn1_-sn11)/(n11i))*((sn_1-sn11)/(n11i+sn-sn1_-sn_1));
					sn11 = n11i;
					return sprob;
				}
				if(n11i==sn11-1){
					sprob *= ((sn11)/(sn1_-n11i))*((sn11+sn-sn1_-sn_1)/(sn_1-n11i));
					sn11 = n11i;
					return sprob;
				}
			}
			sn11 = n11i;
		}
		else{
			sn11 = n11i;
			sn1_=n1_i;
			sn_1=n_1i;
			sn=ni;
		}
		sprob = hyper_323(sn11,sn1_,sn_1,sn);
		return sprob;
	}

	private static double exact(double n11,double n1_,double n_1,double n){
		double p,i,j,prob;
		double max=n1_;
		if(n_1<max) max=n_1;
		double min = n1_+n_1-n;
		if(min<0) min=0;
		if(min==max)
		{
			sless = 1;
			sright= 1;
			sleft = 1;
			slarg = 1;
			return 1;
		}
		prob=hyper0(n11,n1_,n_1,n);
		sleft=0;
		p=hyper(min);
		for(i=min+1; p<0.99999999*prob; i++)
		{
			sleft += p;
			p=hyper(i);
		}
		i--;
		if(p<1.00000001*prob) sleft += p;
		else i--;
		sright=0;
		p=hyper(max);
		for(j=max-1; p<0.99999999*prob; j--)
		{
			sright += p;
			p=hyper(j);
		}
		j++;
		if(p<1.00000001*prob) sright += p;
		else j++;
		if(Math.abs(i-n11)<Math.abs(j-n11)) 
		{
			sless = sleft;
			slarg = 1 - sleft + prob;
		} 
		else 
		{
			sless = 1 - sright + prob;
			slarg = sright;
		}
		return prob;
	}

	private static double exact22(int n11,int n12,int n21,int n22, int lr2t){
		if(n11<0) n11 *= -1;
		if(n12<0) n12 *= -1;
		if(n21<0) n21 *= -1;
		if(n22<0) n22 *= -1; 
		if(n11old==n11 && n12old==n12 && n21old==n21 && n22old==n22) return 1;
		n11old=n11;
		n12old=n12;
		n21old=n21;
		n22old=n22;
		int n1_ = n11+n12;
		int n_1 = n11+n21;
		int n   = n11 +n12 +n21 +n22;
		@SuppressWarnings("unused")
		double prob=exact(n11,n1_,n_1,n);
		left    = sless;
		right   = slarg;
		twotail = sleft+sright;
		if(twotail>1) twotail=1;
		switch (lr2t){
		case 0:
			return left;
		case 1:
			return right;
		case 2:
			return twotail;
		}
		return twotail;
	}
	
	public static double exact22(int n11,int n12,int n21,int n22){
		double d = exact22(n11, n12, n21, n22, 2);
		if (DEBUG) Tools.print("\n<" + n11 + "," + n12 + "," + n21 + "," + n22 + ">");
		return d;//Math.round(d*Math.pow(10, Settings.DEC_PLACES_PVAL));
	}
	
	public static double exact22(int[] n){
		return exact22(n[0],n[1],n[2],n[3]);
	}
	
	public static double exact22total(int[] n){
		return exact22(n[0],n[1]-n[0],n[2],n[3]-n[2]);
	}
	
	public static double exact22total(Double n11,int n12,Double n21,int n22){
		return exact22(n11.intValue(),n12-n11.intValue(),n21.intValue(),n22-n21.intValue());
	}
	
	public static double exact22total(int n11,int n12,Double n21,int n22){
		return exact22(n11,n12-n11,n21.intValue(),n22-n21.intValue());
	}
	
	public static double exact22total(Double n11,int n12,int n21,int n22){
		return exact22(n11.intValue(),n12-n11.intValue(),n21,n22-n21);
	}
	
	public static double exact22total(Double n11,Double n12,Double n21,Double n22){
		return exact22(n11.intValue(),n12.intValue()-n11.intValue(),n21.intValue(),n22.intValue()-n21.intValue());
	}
	
	public static double exact22total(int n11,int n12,int n21,int n22){
		return exact22(n11,n12-n11,n21,n22-n21);
	}
}