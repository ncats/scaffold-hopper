// $Id: Cache.java 3996 2010-02-01 20:19:07Z nguyenda $

package gov.nih.ncgc.util;

import java.util.Collections;
import java.util.Collection;
import java.util.List;
import java.util.ArrayList;
import java.util.Set;
import java.util.concurrent.*;
import java.util.logging.Logger;
import java.util.logging.Level;
import java.lang.ref.SoftReference;
import java.lang.ref.ReferenceQueue;

/**
 * Poor's man caching 
 */
public class Cache<K,V> implements Runnable {
    private static final Logger logger = Logger.getLogger
	(Cache.class.getName());

    final private ConcurrentMap<K, SoftReference<V>> map = 
	new ConcurrentHashMap<K, SoftReference<V>>();

    // reverse lookup
    final private ConcurrentMap<SoftReference<V>, K> rev = 
	new ConcurrentHashMap<SoftReference<V>, K>();

    final private ReferenceQueue<V> queue = new ReferenceQueue<V>();

    public Cache () {
	this (60); // 1min
    }

    public Cache (long refreshRate) {
	ScheduledExecutorService thread = Executors.newScheduledThreadPool(1);
	thread.scheduleAtFixedRate
	    (this, 5, refreshRate, TimeUnit.SECONDS);
    }

    public void run () {
	cleanStaleEntries ();
    }

    public void clear () { 
	map.clear();
	rev.clear();
    }
    public V put (K key, V value) {
	V val = null;
	if (key != null) {
	    //cleanStaleEntries ();
	    SoftReference<V> newRef = new SoftReference<V>(value, queue);
	    SoftReference<V> oldRef = map.put(key, newRef);
	    if (oldRef != null) {
		val = oldRef.get();
		map.remove(rev.remove(oldRef));
	    }
	    rev.put(newRef, key);
	}
	else {
	    val = value;
	}
	return val;
    }
    public V get (K key) {
	if (key != null) {
	    //cleanStaleEntries ();
	    SoftReference<V> ref = map.get(key);
	    if (ref != null) {
		V val = ref.get();
		if (val == null) { // remove this key if its value has been gc
		    map.remove(key);
		    rev.remove(ref);
		}
		return val;
	    }
	}
	return null;
    }

    public V remove (Object key) {
	if (key != null) {
	    //cleanStaleEntries ();
	    SoftReference<V> ref = map.remove(key);
	    if (ref != null) {
		rev.remove(ref);
		V val = ref.get();
		ref.clear();
		return val;
	    }
	}
	return null;
    }
    public int size () { 
	//cleanStaleEntries ();
	return map.size(); 
    }
    public Set<K> keySet () { 
	//cleanStaleEntries ();
	return Collections.unmodifiableSet(map.keySet()); 
    }
    public Collection<V> values () {
	//cleanStaleEntries ();
	List<V> v = new ArrayList<V>();
	for (SoftReference<V> e : map.values()) {
	    v.add(e.get());
	}
	return v;
    }
    public boolean isEmpty () { 
	//cleanStaleEntries ();
	return map.isEmpty(); 
    }
    public boolean containsKey (K key) {
	//cleanStaleEntries ();
	return map.containsKey(key);
    }

    private void cleanStaleEntries () {
	for (SoftReference<V> ref; 
	     (ref = (SoftReference)queue.poll()) != null; ) {
	    K key = rev.remove(ref);
	    if (key != null) {
		map.remove(key);
	    }
	}
    }
}
