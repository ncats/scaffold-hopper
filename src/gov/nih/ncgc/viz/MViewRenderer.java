package gov.nih.ncgc.viz;

import java.awt.*;
import java.awt.geom.AffineTransform;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.awt.geom.RectangularShape;
import java.awt.geom.RoundRectangle2D;
import java.awt.image.BufferedImage;

import chemaxon.marvin.beans.MViewPane;
import chemaxon.marvin.MolPrinter;
import chemaxon.struc.Molecule;
import chemaxon.struc.MolAtom;
import chemaxon.struc.MDocument;
import chemaxon.util.MolHandler;
import chemaxon.formats.MolFormatException;
import chemaxon.marvin.paint.DispOptConsts;

import prefuse.render.LabelRenderer;
import prefuse.render.AbstractShapeRenderer;
import prefuse.visual.VisualItem;
import prefuse.util.ColorLib;
import prefuse.util.FontLib;
import prefuse.util.GraphicsLib;
import prefuse.util.StringLib;
import prefuse.Constants;


public class MViewRenderer extends AbstractShapeRenderer {
    final static int MOL_W = 120;
    final static int MOL_H = 120;
    protected String m_delim = "\n";

    private MolPrinter m_molpainter;
    private String m_molfield = "_structure";
    private String m_textfield = "label";

    private Molecule m_structure = null;
    protected String    m_text; // label text
    private boolean m_molvisible = true;
    private boolean m_textvisible = true;

    protected int m_xAlign = Constants.CENTER;
    protected int m_yAlign = Constants.CENTER;
    protected int m_hTextAlign = Constants.CENTER;
    protected int m_vTextAlign = Constants.CENTER;
    protected int m_hMolAlign = Constants.CENTER;
    protected int m_vMolAlign = Constants.CENTER;
    protected int m_molPos = Constants.TOP;
    
    protected int m_horizBorder = 2;
    protected int m_vertBorder  = 0;
    protected int m_molMargin = 2;
    protected int m_arcWidth    = 0;
    protected int m_arcHeight   = 0;

    protected int m_maxTextWidth = -1;
    
    /** Transform used to scale and position images */
    AffineTransform m_transform = new AffineTransform();

    /** The holder for the currently computed bounding box */
    protected Shape m_bbox  = new Rectangle2D.Double();
    protected Point2D m_pt = new Point2D.Double(); // temp point
    protected Font    m_font; // temp font holder
    protected Dimension m_textDim = new Dimension(); // text width / height
    protected Dimension m_strucDim = new Dimension (MOL_W, MOL_H);
    protected Rectangle m_strucbox = new Rectangle ();
    protected BasicStroke m_stroke = 
	new BasicStroke (1.5f, BasicStroke.CAP_ROUND,BasicStroke.JOIN_ROUND);

    public MViewRenderer () {
	init ();
    }

    public MViewRenderer (String molfield, String textfield) {
	setMolField (molfield);
	setTextField (textfield);
	init ();
    }

    protected void init () {
	m_molpainter = new MolPrinter ();
	m_molpainter.setDispopts(/*DispOptConsts.ATNUM_FLAG | */
				 DispOptConsts.CPK_SCHEME |
				 DispOptConsts.RGROUPS_FLAG |
				 DispOptConsts.QUALITY_MASK);
	m_molpainter.setBackgroundColor(new Color (255, 255, 255, 0));

	//m_molpainter.setColorScheme(DispOptConsts.MONO_SCHEME_S);
    }

    protected BasicStroke getStroke (VisualItem item) {
	return m_stroke;
    }

    public void setShape (Shape shape) { 
        m_bbox = shape; 
        m_strucDim.width = shape.getBounds().width;
        m_strucDim.height = shape.getBounds().height;
    }
    public Shape getShape () { return m_bbox; }

    public void setStroke (BasicStroke stroke) {
	m_stroke = stroke;
    }

    public void setMolField (String label) {
	m_molfield = label;
    }
    public String getMolField () { 
	return m_molfield;
    }
    public void setTextField (String label) {
	m_textfield = label;
    }
    public void setMolVisible (boolean visible) {
	m_molvisible = visible;
    }
    public boolean isMolVisible () { return m_molvisible; }
    public void setTextVisible (boolean visible) {
	m_textvisible = visible;
    }
    public boolean isTextVisible () { return m_textvisible; }

    protected Molecule getMol (VisualItem item) {
	Molecule m = null;

        if (item.canGet(m_molfield, Molecule.class)) {
	    m = (Molecule)item.get(m_molfield);
	}
	else {
	    //item.getTable().addColumn("_structure", Molecule.class);
	    return null;
	}

	if (m == null && item.canGetString(m_molfield) ) {
	    String mol = item.getString(m_molfield);
	    if (mol != null) {
		try {
		    MolHandler mh = new MolHandler (mol);
		    m = mh.getMolecule();
                    if (m.getDim() < 2) 
                        m.clean(2, null);
		    m.dearomatize();
		    item.set(m_molfield, m);
		}
		catch (MolFormatException ex) {
                    ex.printStackTrace();
		}
	    }
	}
        return m;
    }

    private String computeTextDimensions
	(VisualItem item, String text, double size) {

        // put item font in temp member variable
        m_font = item.getFont();
        // scale the font as needed
        if ( size != 1 ) {
            m_font = FontLib.getFont(m_font.getName(), m_font.getStyle(),
                                     size*m_font.getSize());
        }
        
        FontMetrics fm = DEFAULT_GRAPHICS.getFontMetrics(m_font);
        StringBuffer str = null;
        
        // compute the number of lines and the maximum width
        int nlines = 1, w = 0, start = 0, end = text.indexOf(m_delim);
        m_textDim.width = 0;
        String line;
        for ( ; end >= 0; ++nlines ) {
            w = fm.stringWidth(line=text.substring(start,end));
            // abbreviate line as needed
            if ( m_maxTextWidth > -1 && w > m_maxTextWidth ) {
                if ( str == null )
                    str = new StringBuffer(text.substring(0,start));
                str.append(StringLib.abbreviate(line, fm, m_maxTextWidth));
                str.append(m_delim);
                w = m_maxTextWidth;
            } else if ( str != null ) {
                str.append(line).append(m_delim);
            }
            // update maximum width and substring indices
            m_textDim.width = Math.max(m_textDim.width, w);
            start = end+1;
            end = text.indexOf(m_delim, start);
        }
        w = fm.stringWidth(line=text.substring(start));
        // abbreviate line as needed
        if ( m_maxTextWidth > -1 && w > m_maxTextWidth ) {
            if ( str == null )
                str = new StringBuffer(text.substring(0,start));
            str.append(StringLib.abbreviate(line, fm, m_maxTextWidth));
            w = m_maxTextWidth;
        } else if ( str != null ) {
            str.append(line);
        }
        // update maximum width
        m_textDim.width = Math.max(m_textDim.width, w);
        
        // compute the text height
        m_textDim.height = fm.getHeight() * nlines;
        
        return str==null ? text : str.toString();
    }

    protected Shape getRawShape(VisualItem item) {
        m_structure = getMol (item);
	m_text = item.canGetString(m_textfield) 
            ? item.getString(m_textfield) : null;

        double size = item.getSize();
        
        // get image dimensions
        double iw=m_strucDim.getWidth(), ih=m_strucDim.getHeight();
	if (m_structure == null) {
	    iw = ih = 0;
	}

        // get text dimensions
        int tw=0, th=0;
        if ( m_text != null ) {
            m_text = computeTextDimensions(item, m_text, size);
            th = m_textDim.height;
            tw = m_textDim.width;   
        }
        
        // get bounding box dimensions
        double w=0, h=0;
        switch ( m_molPos ) {
        case Constants.LEFT:
        case Constants.RIGHT:
            w = tw + size*(iw +2*m_horizBorder
                   + (tw>0 && iw>0 ? m_molMargin : 0));
            h = Math.max(th, size*ih) + size*2*m_vertBorder;
            break;
        case Constants.TOP:
        case Constants.BOTTOM:
            w = Math.max(tw, size*iw) + size*2*m_horizBorder;
            h = th + size*(ih + 2*m_vertBorder
                   + (th>0 && ih>0 ? m_molMargin : 0));
            break;
        default:
            throw new IllegalStateException(
                "Unrecognized image alignment setting.");
        }

	w = Math.max(5, w);
	h = Math.max(5, h);
        
        // get the top-left point, using the current alignment settings
        getAlignedPoint(m_pt, item, w, h, m_xAlign, m_yAlign);
        
        if ( m_bbox instanceof RoundRectangle2D ) {
            RoundRectangle2D rr = (RoundRectangle2D)m_bbox;
            rr.setRoundRect(m_pt.getX(), m_pt.getY(), w, h,
                            size*m_arcWidth, size*m_arcHeight);
        } 
        else if (m_bbox instanceof RectangularShape) {
            ((RectangularShape)m_bbox).setFrame
                (m_pt.getX(), m_pt.getY(), w, h);
        }
        return m_bbox;
    }
    
    /**
     * Helper method, which calculates the top-left co-ordinate of an item
     * given the item's alignment.
     */
    protected static void getAlignedPoint(Point2D p, VisualItem item, 
					  double w, double h, int xAlign, 
					  int yAlign)
    {
        double x = item.getX(), y = item.getY();
        if ( Double.isNaN(x) || Double.isInfinite(x) )
            x = 0; // safety check
        if ( Double.isNaN(y) || Double.isInfinite(y) )
            y = 0; // safety check
        
        if ( xAlign == Constants.CENTER ) {
            x = x-(w/2);
        } else if ( xAlign == Constants.RIGHT ) {
            x = x-w;
        }
        if ( yAlign == Constants.CENTER ) {
            y = y-(h/2);
        } else if ( yAlign == Constants.BOTTOM ) {
            y = y-h;
        }
        p.setLocation(x,y);
    }

    public void render (Graphics2D g, VisualItem item) {
        RectangularShape shape = (RectangularShape)getShape(item);
        if ( shape == null ) return;
        
        // fill the shape, if requested
        int type = getRenderType(item);
        if ( type==RENDER_TYPE_FILL || type==RENDER_TYPE_DRAW_AND_FILL )
            GraphicsLib.paint
		(g, item, shape, getStroke(item), RENDER_TYPE_FILL);

	
        // now render the image and text
	if (m_text == null && m_structure == null)
	    return;

        double size = item.getSize();
        boolean useInt = 1.5 > Math.max(g.getTransform().getScaleX(),
                                        g.getTransform().getScaleY());
        double x = shape.getMinX() + size*m_horizBorder;
        double y = shape.getMinY() + size*m_vertBorder;


        // render molecule
        if (m_structure != null ) {            
            double w = size * m_strucDim.width;
            double h = size * m_strucDim.height;
            double ix=x, iy=y;

            // determine one co-ordinate based on the image position
            switch ( m_molPos ) {
            case Constants.LEFT:
                x += w + size*m_molMargin;
                break;
            case Constants.RIGHT:
                ix = shape.getMaxX() - size*m_horizBorder - w;
                break;
            case Constants.TOP:
                y += h + size*m_molMargin;
                break;
            case Constants.BOTTOM:
                iy = shape.getMaxY() - size*m_vertBorder - h;
                break;
            default:
                throw new IllegalStateException
		    ("Unrecognized image alignment setting.");
            }
            
            // determine the other coordinate based on image alignment
            switch ( m_molPos ) {
            case Constants.LEFT:
            case Constants.RIGHT:
                // need to set image y-coordinate
                switch ( m_vMolAlign ) {
                case Constants.TOP:
                    break;
                case Constants.BOTTOM:
                    iy = shape.getMaxY() - size*m_vertBorder - h;
                    break;
                case Constants.CENTER:
                    iy = shape.getCenterY() - h/2;
                    break;
                }
                break;
            case Constants.TOP:
            case Constants.BOTTOM:
                // need to set image x-coordinate
                switch ( m_hMolAlign ) {
                case Constants.LEFT:
                    break;
                case Constants.RIGHT:
                    ix = shape.getMaxX() - size*m_horizBorder - w;
                    break;
                case Constants.CENTER:
                    ix = shape.getCenterX() - w/2;
                    break;
                }
                break;
            }

	    AffineTransform tx = g.getTransform();
	    Dimension dim = new Dimension ((int)w, (int)h);

	    g.setRenderingHint(RenderingHints.KEY_ANTIALIASING, 
			       RenderingHints.VALUE_ANTIALIAS_ON);
	    g.setRenderingHint(RenderingHints.KEY_RENDERING, 
				RenderingHints.VALUE_RENDER_QUALITY);

	    m_transform.setTransform(size,0,0,size,ix,iy);
	    AffineTransform save = (AffineTransform)tx.clone();
	    tx.concatenate(m_transform);
	    g.setTransform(tx);

            if (m_structure.hasAtomSet() || m_structure.hasBondSet()) {
                for (MolAtom a : m_structure.getAtomArray()) {
                    int map = a.getAtomMap();
                    if (map > 0) {
                        a.setAtomMap(0);
                        a.setSetSeq(3); // blue
                    }
                }
                m_molpainter.setColorScheme(DispOptConsts.MONO_SCHEME_S);
            }
            else {
                m_molpainter.setColorScheme(DispOptConsts.CPK_SCHEME_S);
            }
            MDocument doc = new MDocument (m_structure);
            doc.setBondSetThickness(0, .125);
            m_molpainter.setDoc(doc);

	    m_strucbox.width = dim.width;
	    m_strucbox.height = dim.height;
	    m_molpainter.setScale(m_molpainter.maxScale(dim));

	    m_molpainter.paint(g, m_strucbox);
	    g.setTransform(save);

            // MolPrinter class turns off antialiasing hints
	    g.setRenderingHint(RenderingHints.KEY_ANTIALIASING, 
			       RenderingHints.VALUE_ANTIALIAS_ON);
	    g.setRenderingHint(RenderingHints.KEY_RENDERING, 
                               RenderingHints.VALUE_RENDER_QUALITY);
        }
        
        // render text
        int textColor = item.getTextColor();
        if ( m_text != null && ColorLib.alpha(textColor) > 0 ) {
            g.setPaint(ColorLib.getColor(textColor));
            g.setFont(m_font);
            FontMetrics fm = DEFAULT_GRAPHICS.getFontMetrics(m_font);

            // compute available width
            double tw;
            switch ( m_molPos ) {
            case Constants.TOP:
            case Constants.BOTTOM:
                tw = shape.getWidth() - 2*size*m_horizBorder;
                break;
            default:
                tw = m_textDim.width;
            }
            
            // compute available height
            double th;
            switch ( m_molPos ) {
            case Constants.LEFT:
            case Constants.RIGHT:
                th = shape.getHeight() - 2*size*m_vertBorder;
                break;
            default:
                th = m_textDim.height;
            }
            
            // compute starting y-coordinate
            y += fm.getAscent();
            switch ( m_vTextAlign ) {
            case Constants.TOP:
                break;
            case Constants.BOTTOM:
                y += th - m_textDim.height;
                break;
            case Constants.CENTER:
                y += (th - m_textDim.height)/2;
            }
            
            // render each line of text
            int lh = fm.getHeight(); // the line height
            int start = 0, end = m_text.indexOf(m_delim);
            for ( ; end >= 0; y += lh ) {
                drawString(g, fm, m_text.substring(start, end), 
			   useInt, x, y, tw);
                start = end+1;
                end = m_text.indexOf(m_delim, start);   
            }
            drawString(g, fm, m_text.substring(start), useInt, x, y, tw);
        }
    
        // draw border
        if (item.isHighlighted() || type==RENDER_TYPE_DRAW || type==RENDER_TYPE_DRAW_AND_FILL) {
            GraphicsLib.paint(g,item,shape,getStroke(item),RENDER_TYPE_DRAW);
        }
    }

    private final void drawString(Graphics2D g, FontMetrics fm, String text,
            boolean useInt, double x, double y, double w)
    {
        // compute the x-coordinate
        double tx;
        switch ( m_hTextAlign ) {
        case Constants.LEFT:
            tx = x;
            break;
        case Constants.RIGHT:
            tx = x + w - fm.stringWidth(text);
            break;
        case Constants.CENTER:
            tx = x + (w - fm.stringWidth(text)) / 2;
            break;
        default:
            throw new IllegalStateException(
                    "Unrecognized text alignment setting.");
        }
        // use integer precision unless zoomed-in
        // results in more stable drawing
        if ( useInt ) {
            g.drawString(text, (int)tx, (int)y);
        } else {
            g.drawString(text, (float)tx, (float)y);
        }
    }

    /**
     * Get the horizontal text alignment within the layout. One of
     * {@link prefuse.Constants#LEFT}, {@link prefuse.Constants#RIGHT}, or
     * {@link prefuse.Constants#CENTER}. The default is centered text.
     * @return the horizontal text alignment
     */
    public int getHorizontalTextAlignment() {
        return m_hTextAlign;
    }
    
    /**
     * Set the horizontal text alignment within the layout. One of
     * {@link prefuse.Constants#LEFT}, {@link prefuse.Constants#RIGHT}, or
     * {@link prefuse.Constants#CENTER}. The default is centered text.
     * @param halign the desired horizontal text alignment
     */
    public void setHorizontalTextAlignment(int halign) {
        if ( halign != Constants.LEFT &&
             halign != Constants.RIGHT &&
             halign != Constants.CENTER )
           throw new IllegalArgumentException(
                   "Illegal horizontal text alignment value.");
        m_hTextAlign = halign;
    }
    
    /**
     * Get the vertical text alignment within the layout. One of
     * {@link prefuse.Constants#TOP}, {@link prefuse.Constants#BOTTOM}, or
     * {@link prefuse.Constants#CENTER}. The default is centered text.
     * @return the vertical text alignment
     */
    public int getVerticalTextAlignment() {
        return m_vTextAlign;
    }
    
    /**
     * Set the vertical text alignment within the layout. One of
     * {@link prefuse.Constants#TOP}, {@link prefuse.Constants#BOTTOM}, or
     * {@link prefuse.Constants#CENTER}. The default is centered text.
     * @param valign the desired vertical text alignment
     */
    public void setVerticalTextAlignment(int valign) {
        if ( valign != Constants.TOP &&
             valign != Constants.BOTTOM &&
             valign != Constants.CENTER )
            throw new IllegalArgumentException(
                    "Illegal vertical text alignment value.");
        m_vTextAlign = valign;
    }
    
    public void setRoundedCorner(int arcWidth, int arcHeight) {
        if ( (arcWidth == 0 || arcHeight == 0) && 
            !(m_bbox instanceof Rectangle2D) ) {
            m_bbox = new Rectangle2D.Double();
        } else {
            if ( !(m_bbox instanceof RoundRectangle2D) )
                m_bbox = new RoundRectangle2D.Double();
            ((RoundRectangle2D)m_bbox)
                .setRoundRect(0,0,10,10,arcWidth,arcHeight);
            m_arcWidth = arcWidth;
            m_arcHeight = arcHeight;
        }
    }
    
    /**
     * Get the image position, determining where the image is placed with
     * respect to the text. One of {@link Constants#LEFT},
     * {@link Constants#RIGHT}, {@link Constants#TOP}, or
     * {@link Constants#BOTTOM}.  The default is left.
     * @return the image position
     */
    public int getMolPosition() {
        return m_molPos;
    }
    
    /**
     * Set the image position, determining where the image is placed with
     * respect to the text. One of {@link Constants#LEFT},
     * {@link Constants#RIGHT}, {@link Constants#TOP}, or
     * {@link Constants#BOTTOM}.  The default is left.
     * @param pos the desired image position
     */
    public void setMolPosition(int pos) {
        if ( pos != Constants.TOP &&
             pos != Constants.BOTTOM &&
             pos != Constants.LEFT &&
             pos != Constants.RIGHT &&
             pos != Constants.CENTER )
           throw new IllegalArgumentException(
                   "Illegal mol position value.");
        m_molPos = pos;
    }
    
    // ------------------------------------------------------------------------
    
    /**
     * Get the horizontal alignment of this node with respect to its
     * x, y coordinates.
     * @return the horizontal alignment, one of
     * {@link prefuse.Constants#LEFT}, {@link prefuse.Constants#RIGHT}, or
     * {@link prefuse.Constants#CENTER}.
     */
    public int getHorizontalAlignment() {
        return m_xAlign;
    }
    
    /**
     * Get the vertical alignment of this node with respect to its
     * x, y coordinates.
     * @return the vertical alignment, one of
     * {@link prefuse.Constants#TOP}, {@link prefuse.Constants#BOTTOM}, or
     * {@link prefuse.Constants#CENTER}.
     */
    public int getVerticalAlignment() {
        return m_yAlign;
    }
    
    /**
     * Set the horizontal alignment of this node with respect to its
     * x, y coordinates.
     * @param align the horizontal alignment, one of
     * {@link prefuse.Constants#LEFT}, {@link prefuse.Constants#RIGHT}, or
     * {@link prefuse.Constants#CENTER}.
     */ 
    public void setHorizontalAlignment(int align) {
        m_xAlign = align;
    }
    
    /**
     * Set the vertical alignment of this node with respect to its
     * x, y coordinates.
     * @param align the vertical alignment, one of
     * {@link prefuse.Constants#TOP}, {@link prefuse.Constants#BOTTOM}, or
     * {@link prefuse.Constants#CENTER}.
     */ 
    public void setVerticalAlignment(int align) {
        m_yAlign = align;
    }
    
    /**
     * Returns the amount of padding in pixels between the content 
     * and the border of this item along the horizontal dimension.
     * @return the horizontal padding
     */
    public int getHorizontalPadding() {
        return m_horizBorder;
    }
    
    /**
     * Sets the amount of padding in pixels between the content 
     * and the border of this item along the horizontal dimension.
     * @param xpad the horizontal padding to set
     */
    public void setHorizontalPadding(int xpad) {
        m_horizBorder = xpad;
    }
    
    /**
     * Returns the amount of padding in pixels between the content 
     * and the border of this item along the vertical dimension.
     * @return the vertical padding
     */
    public int getVerticalPadding() {
        return m_vertBorder;
    }
    
    /**
     * Sets the amount of padding in pixels between the content 
     * and the border of this item along the vertical dimension.
     * @param ypad the vertical padding
     */
    public void setVerticalPadding(int ypad) {
        m_vertBorder = ypad;
    }
    
    /**
     * Get the padding, in pixels, between an image and text.
     * @return the padding between an image and text
     */
    public int getMolTextPadding() {
        return m_molMargin;
    }
    
    /**
     * Set the padding, in pixels, between an image and text.
     * @param pad the padding to use between an image and text
     */
    public void setMolTextPadding(int pad) {
        m_molMargin = pad;
    }
}
