package gov.nih.ncgc.search;


public class SearchParams implements java.io.Serializable {
    private static final long serialVersionUID = 0x07424ca4cf6c8e9al;

    public enum Type {
        Substructure,
            Superstructure,
            Similarity,
            Exact
            };


    // hard limit on screen results
    protected static final int DEFAULT_MATCH_LIMIT = 10000;
    // minimum fingerprint density required to proceed
    protected static final double DEFAULT_MIN_DENSITY = 0.005;
    // milliseconds
    protected static final long DEFAULT_SEARCH_TIMEOUT = 1000*60;
    protected static final double DEFAULT_SIMILARITY = 0.7;

    protected int matchLimit = DEFAULT_MATCH_LIMIT;
    protected double minDensity = DEFAULT_MIN_DENSITY;
    protected long searchTimeout = DEFAULT_SEARCH_TIMEOUT;
    protected double similarity = DEFAULT_SIMILARITY;
    protected Type type;
    protected String rankBy=null;
    protected boolean aromatize = true;
    protected boolean tautomerize = false;

    public SearchParams () {
        this (Type.Substructure);
    }
    public SearchParams (Type type) {
    	this.type=type;
    }
    //TP 08.22.2012: Added specifiable "RANKBY" property for dynamic sorting/paging
    public SearchParams (Type type,String rankBy) {
        this.type = type;
        this.rankBy = rankBy;
    }

    public SearchParams setAromatize (boolean aromatize) { 
        this.aromatize = aromatize;
        return this;
    }
    public boolean getAromatize () { return aromatize; }

    public SearchParams setTautomerize (boolean tautomerize) {
        this.tautomerize = tautomerize;
        return this;
    }
    public boolean getTautomerize () { return tautomerize; }
    
    //TP 08.22.2012: Added setter/getter for rankBy
    public void setRankBy(String rankBy){
    	this.rankBy=rankBy;
    }
    public String getRankBy(){
    	return this.rankBy;
    }

    public static SearchParams substructure () { 
        SearchParams params = new SearchParams (Type.Substructure); 
        return params.setTautomerize(true);
    }
    public static SearchParams superstructure () {
        SearchParams params = new SearchParams (Type.Superstructure);
        return params.setTautomerize(true);
    }
    public static SearchParams similarity () {
        return new SearchParams (Type.Similarity);
    }
    public static SearchParams similarity (double similarity) {
        return new SearchParams (Type.Similarity).setSimilarity(similarity);
    }
    public static SearchParams exact () {
        return new SearchParams (Type.Exact);
    }

    public SearchParams setType (Type type) { 
        this.type = type;
        return this;
    }
    public Type getType () { return type; }

    public SearchParams setMatchLimit (int limit) {
	this.matchLimit = limit;
        return this;
    }
    public int getMatchLimit () { return matchLimit; }

    public SearchParams setMinDensity (double density) {
	this.minDensity = density;
        return this;
    }
    public double getMinDensity () { return minDensity; }

    public SearchParams setTimeout (long timeout) { // in miliseconds
        searchTimeout = timeout;
        return this;
    }
    public long getTimeout () { return searchTimeout; }

    public SearchParams setSimilarity (double similarity) {
        this.similarity = similarity;
        return this;
    }
    public double getSimilarity () { return similarity; }

    public String toString () {
        return "{type="+type+",limit="+matchLimit+",density="+minDensity
            +",timeout="+searchTimeout+",similarity="+similarity+"}";
    }
}
