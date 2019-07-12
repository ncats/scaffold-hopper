// $Id: Base52.java 2950 2009-07-19 07:33:07Z nguyenda $

package gov.nih.ncgc.util;

import java.io.*;

public class Base52 {
    private static final String _A 
	= "AaBbCcDdEeFfGgHhIiJjKkLlMmNnOoPpQqRrSsTtUuVvWwXxYyZz";
    /*
    private static final String _A 
	= "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
    */
    // 2^17 (131072) =~ 52^3 (140608)
    private static String[] TRIGRAM;
    static {
	long start = System.currentTimeMillis();

	// now generate the lut
	TRIGRAM = new String[135252];

	int n = _A.length(), l = 0;
	for (int i = 0; i < n; ++i) {
	    char c = _A.charAt(i);
	    for (int j = 0; j < n; ++j) {
		char b = _A.charAt(j);
		if (c != b) {
		    for (int k = 0; k < n; ++k) {
			char a = _A.charAt(k);
			if (a != b) {
			    TRIGRAM[l++] = String.valueOf(a)
				+ String.valueOf(b)
				+ String.valueOf(c);
			}
		    }
		}
	    }
	}
	/*
	  System.err.println("** Generated base52 LUT in "
	  + String.format("%1$.3f(s)", 
	  (System.currentTimeMillis() 
	  - start)/1000.));
	*/
    }

    public static String encode17 (byte... data) {
	if (data.length < 3) {
	    throw new IllegalArgumentException ("Not enough data specified");
	}
	return encode17(new StringBuffer (3), data).toString();
    }

    public static StringBuffer encode17 (StringBuffer sb, byte... data) {
	int b0 = (data[0]&0xff)|((data[1]&0xff)<<8)|((data[2]&0x01)<<16);
	return sb.append(TRIGRAM[b0]);
    }

    // encode data with 34 bits
    public static String encode34 (byte... data) {
	if (data.length < 5) {
	    throw new IllegalArgumentException ("Not enough data specified");
	}
	return encode34(new StringBuffer (6), data).toString();
    }
    
    public static StringBuffer encode34 (StringBuffer sb, byte... data) {
	int b1 = ((data[2]&0xfe)|((data[3]&0xff)<<8)|((data[4]&0x03)<<16))>>1;
	return encode17(sb, data).append(TRIGRAM[b1]);
    }

    public static String encode51 (byte... data) {
	if (data.length < 7) {
	    throw new IllegalArgumentException ("Not enough data specified");
	}
	return encode51(new StringBuffer (9), data).toString();
    }

    public static StringBuffer encode51 (StringBuffer sb, byte... data) {
	int b2 = ((data[4]&0xfc)|((data[5]&0xff)<<8)|((data[6]&0x07)<<16))>>2;
	return encode34(sb, data).append(TRIGRAM[b2]);
    }

    public static String encode68 (byte... data) {
	if (data.length < 9) {
	    throw new IllegalArgumentException ("Not enough data specified");
	}
	return encode68(new StringBuffer (12), data).toString();
    }

    public static String encode (long x) {
	byte[] b = new byte[9];
	b[8] = 0;
	b[0] = (byte)((x & 0xff00000000000000l) >> 56);
	b[1] = (byte)((x & 0x00ff000000000000l) >> 48);
	b[2] = (byte)((x & 0x0000ff0000000000l) >> 40);
	b[3] = (byte)((x & 0x000000ff00000000l) >> 32);
	b[4] = (byte)((x & 0x00000000ff000000l) >> 24);
	b[5] = (byte)((x & 0x0000000000ff0000l) >> 16);
	b[6] = (byte)((x & 0x000000000000ff00l) >>  8);
	b[7] = (byte)(x & 0x00000000000000ffl);

	return encode68 (b);
    }

    public static StringBuffer encode68 (StringBuffer sb, byte... data) {
	int b3 = ((data[6]&0xf8)|((data[7]&0xff)<<8)|((data[8]&0x0f)<<16))>>3;
	return encode51(sb, data).append(TRIGRAM[b3]);
    }

    public static String encode85 (byte... data) {
	if (data.length < 11) {
	    throw new IllegalArgumentException ("Not enough data specified");
	}
	return encode85(new StringBuffer (15), data).toString();
    }

    public static StringBuffer encode85 (StringBuffer sb, byte... data) {
	int b4=((data[8]&0xf0)|((data[9]&0xff)<<8)|((data[10]&0x1f)<<16))>>4;
	return encode68(sb, data).append(TRIGRAM[b4]);
    }

    public static String encode102 (byte... data) {
	if (data.length < 13) {
	    throw new IllegalArgumentException ("Not enough data specified");
	}
	return encode102(new StringBuffer (18), data).toString();
    }

    public static StringBuffer encode102 (StringBuffer sb, byte... data) {
	int b5 = ((data[10]&0xe0)|((data[11]&0xff)<<8)
		  |((data[12]&0x3f)<<16))>>5;
	return encode85(sb, data).append(TRIGRAM[b5]);
    }

    public static String encode119 (byte... data) {
	if (data.length < 15) {
	    throw new IllegalArgumentException ("Not enough data specified");
	}
	return encode119(new StringBuffer (21), data).toString();
    }

    public static StringBuffer encode119 (StringBuffer sb, byte... data) {
	int b6 = ((data[12]&0xc0)|((data[13]&0xff)<<8)
		  | ((data[14]&0x7f)<<16))>>6;
	return encode102(sb, data).append(TRIGRAM[b6]);
    }

    public static String encode136 (byte... data) {
	if (data.length < 17) {
	    throw new IllegalArgumentException ("Not enough data specified");
	}
	return encode136(new StringBuffer (24), data).toString();
    }

    public static StringBuffer encode136 (StringBuffer sb, byte... data) {
	int b7 = ((data[14]&0x80)|((data[15]&0xff)<<8)
		  | ((data[16]&0xff)<<16))>>7;
	return encode119(sb, data).append(TRIGRAM[b7]);
    }

    public static String encode153 (byte... data) {
	if (data.length < 20) {
	    throw new IllegalArgumentException ("Not enough data specified");
	}
	// 27 = 3*(153/17)
	return encode153(new StringBuffer (27), data).toString();
    }

    public static StringBuffer encode153 (StringBuffer sb, byte... data) {
	int b8 = ((data[17] & 0xff) | ((data[18] & 0xff) << 8) 
		  | ((data[19] & 0x01)<<16));
	return encode136(sb, data).append(TRIGRAM[b8]);
    }

    public static void encode153 (int[] indexes, byte... data) {
	int b0 = (data[0]&0xff)|((data[1]&0xff)<<8)|((data[2]&0x01)<<16);
	int b1 = ((data[2]&0xfe)|((data[3]&0xff)<<8)|((data[4]&0x03)<<16))>>1;
	int b2 = ((data[4]&0xfc)|((data[5]&0xff)<<8)|((data[6]&0x07)<<16))>>2;
	int b3 = ((data[6]&0xf8)|((data[7]&0xff)<<8)|((data[8]&0x0f)<<16))>>3;
	int b4 =((data[8]&0xf0)|((data[9]&0xff)<<8)|((data[10]&0x1f)<<16))>>4;
	int b5 = ((data[10]&0xe0)|((data[11]&0xff)<<8)
		  |((data[12]&0x3f)<<16))>>5;
	int b6 = ((data[12]&0xc0)|((data[13]&0xff)<<8)
		  | ((data[14]&0x7f)<<16))>>6;
	int b7 = ((data[14]&0x80)|((data[15]&0xff)<<8)
		  | ((data[16]&0xff)<<16))>>7;
	int b8 = ((data[17] & 0xff) | ((data[18] & 0xff) << 8) 
		  | ((data[19] & 0x01)<<16));
	indexes[0] = b0;
	indexes[1] = b1;
	indexes[2] = b2;
	indexes[3] = b3;
	indexes[4] = b4;
	indexes[5] = b5;
	indexes[6] = b6;
	indexes[7] = b7;
	indexes[8] = b8;
    }

    static void generateLUT () throws Exception {
	int n = _A.length(), l = 0;
	for (int i = 0; i < n; ++i) {
	    char c = _A.charAt(i);
	    for (int j = 0; j < n; ++j) {
		char b = _A.charAt(j);
		if (c != b) {
		    for (int k = 0; k < n; ++k) {
			char a = _A.charAt(k);
			if (a != b) {
			    ++l;
			}
		    }
		}
	    }
	}

	System.out.println("lut size: " + l);
	String[] lut = new String[l];
	l = 0;
	for (int i = 0; i < n; ++i) {
	    char c = _A.charAt(i);
	    for (int j = 0; j < n; ++j) {
		char b = _A.charAt(j);
		if (c != b) {
		    for (int k = 0; k < n; ++k) {
			char a = _A.charAt(k);
			if (a != b) {
			    lut[l++] = String.valueOf(c)
				+ String.valueOf(b) 
				+ String.valueOf(a);
			}
		    }
		}
	    }
	}
	
	java.io.ObjectOutputStream oos = new java.io.ObjectOutputStream
	    (new java.io.FileOutputStream ("base52_lut.obj"));
	oos.writeObject(lut);
	oos.close();
    }

    public static void main (String[] argv) throws Exception {
	for (long x = 0; x < 999; ++x) {
	    System.out.printf("%1$s %2$3d\n", encode (x), x);
	}
    }

    public static void main0 (String[] argv) throws Exception {
	java.security.SecureRandom rand = new java.security.SecureRandom();
	byte[] b153= new byte[20];
	byte[] b136= new byte[17];
	byte[] b119= new byte[15];
	byte[] b102= new byte[13];
	byte[] b85 = new byte[11];
	byte[] b68 = new byte[9];
	byte[] b51 = new byte[7];
	byte[] b34 = new byte[5];

	byte[] mesg = new byte[20];
	java.security.MessageDigest digest = 
	    java.security.MessageDigest.getInstance("SHA1");
	byte[] hash;

	int[][] distribution = new int[9][135252];
	int[] indexes = new int[9];

	java.util.Set<String> set = new java.util.HashSet<String>();
	for (int i = 0; i < 1000000; ++i) {
	    rand.nextBytes(mesg);
	    /*
	      digest.reset();
	      hash = digest.digest(mesg);
	    */
	    hash = mesg;

	    System.arraycopy(hash, 0, b153, 0, b153.length);
	    System.arraycopy(hash, 0, b136, 0, b136.length);
	    System.arraycopy(hash, 0, b119, 0, b119.length);
	    System.arraycopy(hash, 0, b102, 0, b102.length);
	    System.arraycopy(hash, 0, b85, 0, b85.length);
	    System.arraycopy(hash, 0, b68, 0, b68.length);
	    System.arraycopy(hash, 0, b51, 0, b51.length);
	    System.arraycopy(hash, 0, b34, 0, b34.length);
	    /*
	      System.out.println(encode34(b34) + " " 
	      + encode51(b51) + " " 
	      + encode68(b68) + " " 
	      + encode85(b85) + " "
	      + encode102(b102) + " "
	      + encode119(b119) + " "
	      + encode136(b136) + " "
	      + encode153(b153));
	    */
	    encode153 (indexes, b153);
	    for (int j = 0; j < indexes.length; ++j) {
		++distribution[j][indexes[j]];
	    }
	}

	for (int i = 0; i < indexes.length; ++i) {
	    int count = 0;
	    int minBucket = -1, maxBucket = -1;
	    int minBucketCount = Integer.MAX_VALUE;
	    int maxBucketCount = 0;
	    for (int j = 0; j < distribution[i].length; ++j) {
		int k = distribution[i][j];
		if (k > 0) {
		    //System.out.println(i + " " + distribution[i]);
		    if (k < minBucketCount) {
			minBucketCount = k;
			minBucket = j;
		    }
		    if (k > maxBucketCount) {
			maxBucketCount = k;
			maxBucket = j;
		    }
		    ++count;
		}
	    }
	    System.out.println("** Bit " + i);
	    System.out.println("hash value distributed " + count 
			       + " across " + distribution[i].length 
			       + " possible buckets");
	    System.out.println("Min bucket: " + minBucket + " " 
			       + TRIGRAM[minBucket] + " size=" 
			       + minBucketCount);
	    System.out.println("Max bucket: " + maxBucket + " "
			       + TRIGRAM[maxBucket] + " size=" 
			       + maxBucketCount);
	}
    }
}
