package gov.nih.ncgc.rgroup;

import java.awt.Component;
import java.awt.Color;
import javax.swing.*;
import javax.swing.table.TableCellRenderer;

import chemaxon.struc.Molecule;
import chemaxon.struc.MolAtom;
import chemaxon.marvin.beans.MViewPane;

public class MolCellRenderer  implements TableCellRenderer {
    MViewPane mview = new MViewPane ();
    JPanel blank = new JPanel ();
    
    public MolCellRenderer () {
    }

    public Component getTableCellRendererComponent 
        (JTable table, Object value, boolean selected, 
         boolean focus, int row, int column) {
	    
        Color bg = selected ?table.getSelectionBackground() 
            : table.getBackground();
        mview.setMolbg(bg);
	    
        if (value instanceof Molecule) {
            Molecule m = (Molecule)value;
		
            if (m != null) {
                //TODO: Fix this trick to be more transparent
                if(m.getName().startsWith("<html>")){
                    return new JLabel(m.getName());
                }
                m.dearomatize();
            }
            mview.setM(0, m);
        }
        else if (value != null) {
            mview.setM(0, value.toString());
        }
        else {
            mview.setM(0, (Molecule)null);
            blank.setBackground(bg);
            return blank;
        }

        return mview;
    }
}
