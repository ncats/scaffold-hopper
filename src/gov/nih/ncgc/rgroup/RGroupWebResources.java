
package gov.nih.ncgc.rgroup;

import java.net.URL;
import java.util.logging.Logger;
import java.util.logging.Level;


public class RGroupWebResources {
    static final Logger logger = Logger.getLogger
        (RGroupWebResources.class.getName());

    static final String DEFAULT_HOST = "https://tripod.nih.gov";
    static final String DEFAULT_BASE = "/chembl/v23";

    protected String host;
    protected String base;

    public RGroupWebResources () {
        this (DEFAULT_HOST, DEFAULT_BASE);
    }

    public RGroupWebResources (String host, String base) {
        this.host = host;
        this.base = base;
    }

    public void setHost (String host) { this.host = host; }
    public String getHost () { return host; }
    public void setBase (String base) { this.base = base; }
    public String getBase () { return base; }

    public String getResource () { return host + "/"+base; }
    public String getCompoundResource (String arg) {
        return getResource()+"/compounds/"+arg;
    }
    public URL getCompoundURL (String arg) {
        String url = getCompoundResource (arg);
        try {
            return new URL (url);
        }
        catch (Exception ex) {
            logger.warning("Bogus compound url: "+url);
        }
        return null;
    }
    public String getDocumentResource (String arg) {
        return getResource()+"/documents/"+arg;
    }
    public URL getDocumentURL (String arg) {
        String url = getDocumentResource (arg);
        try {
            return new URL (url);
        }
        catch (Exception ex) {
            logger.warning("Bogus document url: "+url);
        }
        return null;
    }
}
