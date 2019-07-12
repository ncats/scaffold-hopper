package gov.nih.ncgc.rgroup;

import java.awt.Component;
import java.awt.Color;
import javax.swing.*;
import javax.swing.table.*;

import chemaxon.struc.Molecule;
import chemaxon.struc.MolAtom;
import chemaxon.marvin.beans.MViewPane;


public class MolCellEditor 
    extends AbstractCellEditor implements TableCellEditor {
    MViewPane mview = new MViewPane ();

    public MolCellEditor () {
    }

    public Component getTableCellEditorComponent 
        (JTable table, Object value, boolean selected, 
         int row, int column) {
        
        mview.setM(0, (Molecule)value);
        return mview;
    }

    public Object getCellEditorValue () {
        if (mview.getVisibleCellCount() == 0) {
            return null;
        }

        Molecule m = mview.getM(0);
        if (m != null) {
            mview.applyRotationMatrices();
        }
        return m;
    }
}
