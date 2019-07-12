package gov.nih.ncgc.rgroup;

import java.beans.*;
import java.io.*;
import java.net.*;
import javax.net.ssl.*;
import java.util.*;
import java.util.logging.Logger;
import java.util.logging.Level;
import java.util.concurrent.*;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.*;

import javax.swing.text.AttributeSet;
import javax.swing.text.html.HTML;

import org.jdesktop.jxlayer.JXLayer;
import org.jdesktop.jxlayer.plaf.LayerUI;
import org.jdesktop.jxlayer.plaf.AbstractLayerUI;
import org.jdesktop.swingworker.SwingWorker;
import org.jdesktop.swingx.JXTaskPane;
import org.jdesktop.swingx.JXTaskPaneContainer;
import org.jdesktop.swingx.JXErrorPane;

import org.codehaus.jackson.map.ObjectMapper;
import org.codehaus.jackson.node.ArrayNode;
import org.codehaus.jackson.JsonNode;

import chemaxon.struc.Molecule;
import chemaxon.struc.MolBond;
import chemaxon.struc.MolAtom;

import gov.nih.ncgc.viz.LabelRenderer;
import gov.nih.ncgc.search.SearchParams;
import gov.nih.ncgc.util.DummySSLSocketFactory;

public class RGroupDocsPane extends JPanel implements HyperlinkListener {
    static final Logger logger = 
        Logger.getLogger(RGroupDocsPane.class.getName());

    static final String CSS = "<body style=\"font: 10px 'Lucida Grande',"
        +"Verdana, Helvetica, Arial, Geneva, sans-serif; color: #333;\">";

    class DocsLoader extends SwingWorker<Throwable, JsonNode> {
        URL url;
        ArrayList<JsonNode> docs = new ArrayList<JsonNode>();

        DocsLoader () {
        }

        DocsLoader (URL url) {
            this.url = url;
        }

        public void setURL (URL url) { this.url = url; }

        @Override
        protected Throwable doInBackground () {
            setMessage ("Loading references...");
            try {
                logger.info("Loading "+url);

                ObjectMapper mapper = new ObjectMapper ();
                URLConnection con = url.openConnection();
                if ("https".equals(url.getProtocol())) {
                    ((HttpsURLConnection)con).setSSLSocketFactory
                        (new DummySSLSocketFactory ());
                }
                
                JsonNode node = mapper.readTree(con.getInputStream());
                logger.info("## JSON: "+node);
                if (node.isArray()) {
                    for (Iterator<JsonNode> iter = node.getElements(); 
                         iter.hasNext(); ) {
                        if (Thread.currentThread().isInterrupted())
                            return null;

                        JsonNode n = iter.next();
                        docs.add(n);
                        //publish (mapper.readValue(n, Map.class));
                        //publish (n);
                    }
                }
                else {
                    docs.add(node);
                    //publish (mapper.readValue(node, Map.class));
                    //publish (node);
                }
            }
            catch (EOFException ex) {
                // empty stream...
            }
            catch (Exception ex) {
                ex.printStackTrace();
                return ex;
            }
            return null;
        }

        @Override
        protected void process (JsonNode... docs) {
            for (JsonNode d : docs) {
                //logger.info("** doc "+d);
                addDoc (d);
            }
        }

        @Override
        protected void done () {
            logger.info("## "+url+"..."+docs.size()+" documents loaded!");
            try {
                Throwable t = get ();
                if (t != null) {
                    JXErrorPane.showDialog(t);
                }
                else {
                    for (JsonNode d : docs) 
                        addDoc (d);
                    fireDocsPropertyChange ();
                }
            }
            catch (Exception ex) {
                JXErrorPane.showDialog(ex);
            }
            clearMessage ();
        }
    }

    class ScaffoldLoader extends DocsLoader {
        ScaffoldLoader (Molecule scaffold) {
            this (scaffold, 100);
        }

        ScaffoldLoader (Molecule scaffold, int max) {
            try {
                Molecule q = scaffold.cloneMolecule();
                q.aromatize();
                for (MolBond bond : q.getBondArray()) {
                    // turn off Hs
                    if (bond.getType() == MolBond.AROMATIC) {
                        bond.getAtom1().setQProp("H", -1);
                        bond.getAtom2().setQProp("H", -1);
                    }
                }
                setURL (getWS().getCompoundURL
                        (URLEncoder.encode(getQuery (q), "utf-8")
                         +"/documents?max="+Math.max(max, 1)));
            }
            catch (Exception ex) {
                ex.printStackTrace();
            }
        }
    }

    class SearchLoader extends DocsLoader {
        SearchLoader (String query) {
            this (query, 100);
        }

        SearchLoader (String query, int max) {
            try {
                if (query.indexOf("\"") < 0) {
                    query = URLEncoder.encode(query, "utf-8");
                }
                else {
                    query = query.replaceAll("\"", "%22")
                        .replaceAll("\\[", "%5b").replaceAll("\\]", "%5d");
                }
                setURL (getWS().getDocumentURL
                        (query+"/search?max="+Math.max(max, 1)));
            }
            catch (Exception ex) {
                ex.printStackTrace();
            }
        }
    }

    class StructureSearchLoader extends DocsLoader {
        StructureSearchLoader (Molecule mol, SearchParams params) {
            this (mol, 100, params);
        }

        StructureSearchLoader (Molecule mol, int max, SearchParams params) {
            try {
                String arg = ""; // addition arguments
                String type = "";
                switch (params.getType()) {
                case Substructure: type = "sub"; break;
                case Superstructure: type = "super"; break;
                case Exact: type = "exact"; break;
                case Similarity: 
                    type = "sim"; 
                    arg = "&cutoff="+params.getSimilarity();
                    break;
                }
                    
                setURL (getWS().getCompoundURL
                        (//URLEncoder.encode(getQuery (mol), "utf-8")
                         getQuery (mol)
                         +"/"+type+"?max="+Math.max(max,1)
                         +"&format=pubmed"+arg));
            }
            catch (Exception ex) {
                ex.printStackTrace();
            }
        }
    }

    protected RGroupWebResources webservice = new RGroupWebResources ();
    protected JXTaskPaneContainer taskpane;
    protected String message;
    protected PropertyChangeSupport pcs = new PropertyChangeSupport (this);

    protected java.util.List<JsonNode> docs = new ArrayList<JsonNode>();
    protected LinkedList<Collection<JsonNode>> history = 
        new LinkedList<Collection<JsonNode>>();
    protected LinkedList<Collection<JsonNode>> forward = 
        new LinkedList<Collection<JsonNode>>();

    protected Future future;
    protected ExecutorService threadPool = Executors.newSingleThreadExecutor();

    public RGroupDocsPane () {
        super (new BorderLayout ());
        
        taskpane = new JXTaskPaneContainer ();

        JPanel panel = new JPanel (new BorderLayout ());
        panel.add(new JScrollPane (taskpane));

        final LabelRenderer labeler = new LabelRenderer ();
        labeler.setLabelFilled(false);
        JXLayer layer = new JXLayer (panel, new AbstractLayerUI () {
                @Override
                public long getLayerEventMask() { 
                    /*
                     * this prevents the layer from intercepting the input 
                     * masks such as keyboard and mouse events.  we have to
                     * do this because turning on the masks will mess around 
                     * with the AWT's default toolkit event listener that 
                     * prevents it to run within webstart sandbox.
                     */
                    return 0; 
                }

                @Override
                protected void paintLayer (Graphics2D g2, JXLayer layer) {
                    super.paintLayer(g2, layer);
                    String mesg = getMessage ();
                    if (mesg != null) {
                        Rectangle r = layer.getBounds();
                        labeler.renderLabel(g2, r.height/2, 
                                            5, 5, r.width, r.height, 
                                            mesg);
                    }
                }
            });
        add (layer);
    }

    
    static String getQuery (Molecule mol) {
        String smarts = mol.toFormat("smarts");
        StringBuilder query = new StringBuilder ();
        for (int i = 0; i < smarts.length(); ++i) {
            char ch = smarts.charAt(i);
            switch (ch) {
            case '/': case '\\':
                break;
            case '+':
                query.append("%2b");
                break;
            case '#':
                query.append("%23");
                break;
            default:
                query.append(ch);
            }
        }
        return query.toString();
    }

    JXTaskPane createTaskPane (JsonNode doc, boolean collapse) {
        JXTaskPane tp = new JXTaskPane ();
        tp.setLayout(new BorderLayout ());
        JsonNode title = doc.get("title");
        tp.setTitle("["+doc.get("year").asText()+"] "
                    +((title == null || title.isNull()) 
                      ? "(no title)" : title.asText()));;

        tp.setCollapsed(collapse);
        tp.setScrollOnExpand(true);

        StringBuilder sb = new StringBuilder (CSS);
        sb.append("<table>");
        sb.append("<tr valign=\"top\">"
                  +"<td align=\"right\"><b>Title</b></td><td>"
                  +"<a href=\"rgroup://load/"+doc.get("pubmedId")+"\">"
                  +doc.get("title").asText()+"</a></td></tr>");
        JsonNode doi = doc.get("doi");
        if (doi != null && !doi.isNull()) {
            sb.append("<tr valign=\"top\">"
                      +"<td align=\"right\"><b>DOI</b></td><td>"
                      +"<a href=\"http://dx.doi.org/"+doi.asText()
                      +"\">"+doi.asText()+"</a></td></tr>");
        }
        sb.append("<tr valign=\"top\">"
                  +"<td align=\"right\"><b>PubMed</b></td><td>" 
                  +"<a href=\"http://www.ncbi.nlm.nih.gov/pubmed/?term="+doc.get("pubmedId").asLong()+"\">"
                  +doc.get("pubmedId").asLong()+"</a></td></tr>");

        JsonNode cnt = doc.get("compoundCount");
        if (cnt != null && !cnt.isNull()) {
            String link = "<a href=\"rgroup://load/"+doc.get("pubmedId")+"\">";
            sb.append("<tr valign=\"top\">"
                      +"<td align=\"right\">"
                      +"<b>Compounds</b></a></td><td>"
                      +link+cnt.asInt()+"</a></td></tr>");
        }

        JsonNode abs = doc.get("abstractText");
        if (abs != null && !abs.isNull()) {
            sb.append("<tr valign=\"top\">"
                      +"<td align=\"right\"><b>Abstract</b></td><td>" 
                      +abs.asText()+"</td></tr>");
        }

        JsonNode authors = doc.get("authors");
        if (authors != null && !authors.isNull()) {
            StringBuilder html = new StringBuilder ();
            for (int i = 0; i < authors.size(); ++i) {
                if (i > 0) html.append("; ");
                html.append(authors.get(i).asText());
            }
            sb.append("<tr valign=\"top\">"
                      +"<td align=\"right\"><b>Authors</b></td><td>" 
                      +html+"</td></tr>");
        }

        JsonNode affiliation = doc.get("affiliation");
        if (affiliation != null && !affiliation.isNull()) {
            sb.append("<tr valign=\"top\">"
                      +"<td align=\"right\"><b>Affiliation</b></td><td>" 
                      +affiliation.asText()+"</td></tr>");
        }

        JsonNode journal = doc.get("journal");
        if (journal != null && !journal.isNull()) {
            sb.append("<tr valign=\"top\">"
                      +"<td align=\"right\"><b>Journal</b></td><td>" 
                      +journal.asText()+"</td></tr>");            
        }

        JsonNode targets = doc.get("targets");
        if (targets != null && !targets.isNull()) {
            StringBuilder html = new StringBuilder ();
            int size = Math.min(targets.size(), 10);
            for (int i = 0; i < size; ++i) {
                JsonNode t = targets.get(i);
                int tid = t.get("tid").asInt();
                String name = t.get("name").asText();
                String organism = t.get("organism").asText();
                try {
                    html.append("<a href=\"rgroup://tid:"+tid+"/"
                                +URLEncoder.encode(name, "utf8")+"\">"
                                +name+"</a>&nbsp;");
                }
                catch (Exception ex) {
                    ex.printStackTrace();
                }
            }
            sb.append("<tr valign=\"top\">"
                      +"<td align=\"right\"><b>Targets</b></td><td>" 
                      +html+"</td></tr>");
        }

        JsonNode terms = doc.get("meSHTerms");
        if (terms != null && !terms.isNull()) {
            StringBuilder html = new StringBuilder ();
            for (int i = 0; i < terms.size(); ++i) {
                JsonNode n = terms.get(i);
                String text = n.asText();
                int pos = text.indexOf('*');
                if (pos > 0) {
                    text = text.substring(0, pos);
                }
                try {
                    html.append("<a href=\"rgroup://"+(pos>0?"majr":"mh")+"/"
                                +URLEncoder.encode(text, "utf8")
                                +"\">"+text+"</a>&nbsp;");
                }
                catch (Exception ex) {
                    ex.printStackTrace();
                }
            }
            sb.append("<tr valign=\"top\"><td align=\"right\">"
                      +"<b>MeSH Terms</b></td><td>"+html+"</td></tr>");
        }

        sb.append("</table>");

            
        JEditorPane ep = new JEditorPane ();
        ep.setContentType("text/html");
        ep.setEditable(false);
        ep.setOpaque(false);
        ep.addHyperlinkListener(this);
        ep.setText(sb.toString());

        tp.add(ep);

        logger.info("doc: "+doc.get("title"));
        return tp;
    }

    public RGroupWebResources getWS () { return webservice; }
    public void setWS (RGroupWebResources webservice) {
        this.webservice = webservice;
    }

    public void hyperlinkUpdate (HyperlinkEvent e) {
        AttributeSet attrs = (AttributeSet)e.getSourceElement()
            .getAttributes().getAttribute(HTML.Tag.A);
        String href = (String)attrs.getAttribute(HTML.Attribute.HREF);
        try {
            URI uri = new URI (href);
            String action = uri.getHost();
            // skip /
            String id = uri.getPath().substring(1); 
            String scheme = uri.getScheme();

            logger.info("## uri: "+uri);
            if (e.getEventType() == HyperlinkEvent.EventType.ACTIVATED) {
                logger.info("action="+action+" id="+id);

                if (scheme.startsWith("http")) {
                    logger.info("## launching browser url "+uri);
                    try {
                        Desktop.getDesktop().browse(uri);
                    }
                    catch (Exception ex) {
                        logger.log(Level.SEVERE, 
                                   "Can't launch external browser "
                                   +uri, ex);
                    }
                }
                else if (action.equalsIgnoreCase("tid")) {
                    id = String.valueOf(uri.getPort());
                }
                pcs.firePropertyChange
                    (new PropertyChangeEvent (this, action, null, id));
            }
            else if (e.getEventType() == HyperlinkEvent.EventType.ENTERED) {
                JComponent c = (JComponent)e.getSource();
                if (scheme.startsWith("http")) {
                    c.setToolTipText("Launch browser with "+uri);
                }
                else if (action.equalsIgnoreCase("load")) {
                    c.setToolTipText("Generate scaffolds for document "+id);
                }
                else if (action.equalsIgnoreCase("mh")
                         || action.equalsIgnoreCase("majr")) {
                    c.setToolTipText("Load MeSH term "+id);
                }
                else if (action.startsWith("tid")) {
                    c.setToolTipText("Load target \""
                                     +URLDecoder.decode(id, "utf8")+"\"");
                }
                else {
                    c.setToolTipText(null);
                }
            }
        }
        catch (Exception ex) {
            logger.log(Level.SEVERE, "Can't parse href \""+href+"\"", ex);
        }
    }

    protected void addDoc (JsonNode doc) {
        docs.add(doc);
        taskpane.add(createTaskPane (doc, taskpane.getComponentCount() > 0));
        repaint ();
    }

    public void setMessage (String message) {
        this.message = message;
        repaint ();
    }
    public void clearMessage () { 
        setMessage (null); 
    }
    public String getMessage () { return message; }

    protected void load (DocsLoader loader) {
        if (!docs.isEmpty()) {
            history.push(new ArrayList<JsonNode>(docs));
            docs.clear();
            taskpane.removeAll();
        }

        if (future != null) {
            logger.info("** canceling current task "+future);
            future.cancel(true);
        }
        future = threadPool.submit(loader);
    }

    public void search (Molecule scaffold) {
        load (new ScaffoldLoader (scaffold));
    }

    public void search (Molecule mol, SearchParams params) {
        load (new StructureSearchLoader (mol, params));
    }

    public void search (String query) {
        load (new SearchLoader (query));
    }

    public void load (URL url) {
        load (new DocsLoader (url));
    }

    public void clear () {
        taskpane.removeAll(); // remove all components
        docs.clear();
        repaint ();
    }

    public void pop () {
        if (history.isEmpty())
            return;

        taskpane.removeAll();
        docs.clear();

        Collection<JsonNode> dd = history.pop();
        for (JsonNode d : dd) {
            addDoc (d);
        }
        forward.push(dd);
        fireDocsPropertyChange ();
    }

    public void push () {
        if (forward.isEmpty())
            return;

        taskpane.removeAll();
        docs.clear();

        Collection<JsonNode> dd = forward.pop();
        for (JsonNode d : dd) {
            addDoc (d);
        }
        history.push(dd);
        fireDocsPropertyChange ();
    }        

    public Collection<JsonNode> getDocs () { return docs; }

    protected void fireDocsPropertyChange () {
        pcs.firePropertyChange(new PropertyChangeEvent 
                               (this, "docs", null, 
                                docs.toArray(new JsonNode[0])));
    }

    public void addPropertyChangeListener (PropertyChangeListener l) {
        pcs.addPropertyChangeListener(l);
    }
    public void removePropertyChangeListener (PropertyChangeListener l) {
        pcs.removePropertyChangeListener(l);
    }
}
