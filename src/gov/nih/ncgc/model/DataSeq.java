// $Id: DataSeq.java 3464 2009-10-26 18:45:16Z nguyenda $

package gov.nih.ncgc.model;

/**
 * a simple interface for providing data sequence
 */
public interface DataSeq<T> {
    public int size ();
    public T get (int index);
}
