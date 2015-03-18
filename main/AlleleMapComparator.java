package main;

import java.util.Comparator;
import java.util.HashMap;

public class AlleleMapComparator implements Comparator {
	private HashMap <String, Double[]>alleleMap;
	private int comparedFeature;
	int order;
	
	public AlleleMapComparator(HashMap <String, Double[]>alleleMap, int comparedFeature, int order){
		super();
		this.alleleMap = alleleMap;
		this.comparedFeature = comparedFeature;
		this.order = order;
	}
	
	public int compare(Object allele1, Object allele2) {
		Double i1, i2;
		int ret = -1;

		i1 = alleleMap.get(allele1)[comparedFeature];
		i2 = alleleMap.get(allele2)[comparedFeature];
		if (order == 0 )
			ret = i1.compareTo(i2);
		else
			ret = i2.compareTo(i1);

		if (ret == 0){
			ret = ((String)allele1).compareTo((String)allele2);
		}
		return (ret == 0 ? 1 : ret);
	}
}