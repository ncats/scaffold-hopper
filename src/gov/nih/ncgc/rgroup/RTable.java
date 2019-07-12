package gov.nih.ncgc.rgroup;

import java.util.*;

import java.awt.Component;
import java.awt.Color;
import java.awt.BorderLayout;
import javax.swing.*;
import javax.swing.border.*;
import javax.swing.table.*;

import chemaxon.struc.Molecule;

public class RTable extends JTable {
    static public class TableHeaderRenderer 
	extends JPanel implements TableCellRenderer {
	JTextArea area;
	JLabel sortLabel;
	Color color;

	public TableHeaderRenderer () {
	    color = getBackground ();

	    area = new JTextArea ();
	    area.setOpaque(true);
	    area.setRows(2);
	    area.setLineWrap(true);
	    area.setWrapStyleWord(true);
	    area.setBackground(color);

	    sortLabel = new JLabel ();
	    setLayout (new BorderLayout (1,0));
	    add (area);
	    add (sortLabel, BorderLayout.EAST);
	    setBorder (BorderFactory.createBevelBorder(BevelBorder.RAISED));
	}

	public Color getColor () { return color; }
	public void setColor (Color color) { 
	    this.color = color;
	    area.setBackground(color);
	    setBackground (color);
	}

	public Component getTableCellRendererComponent 
	    (JTable table, Object value, boolean isSelected,
	     boolean hasFocus, int row, int column) {

	    area.setText(value == null || column < 0 
			 ? null : value.toString());
	    return this;
	}
    }

    public RTable () {
        this (new DefaultTableModel ());
    }
    public RTable (TableModel model) {
        super (model);
	setDefaultRenderer (Molecule.class, new MolCellRenderer());
	setDefaultEditor (Molecule.class, new MolCellEditor());
	setDefaultRenderer (ScaffoldData.class, 
                            new ScaffoldDataCellRenderer ());
	getTableHeader().setDefaultRenderer
	    (new TableHeaderRenderer ());
	setRowHeight(100);
	setAutoResizeMode(JTable.AUTO_RESIZE_OFF);

    	setAutoCreateRowSorter (true);
    	TableRowSorter trs = new TableRowSorter(model);
    	for(int i=0; i < model.getColumnCount(); i++){
            if(model.getColumnClass(i).isAssignableFrom(Molecule.class)){
                trs.setComparator
                    (i, new Comparator<Molecule>(){
                        @Override
                        public int compare(Molecule arg0, Molecule arg1) {
                            if(arg0==null)return -1;
                            if(arg1==null)return 1;
                            return Double.compare(arg0.getMass(),arg1.getMass());						
                        }});
            }
    	}
    	
        setRowSorter(trs);
    }

    @Override
    public void createDefaultColumnsFromModel () {
        TableColumnModel columns = getColumnModel ();
        TableModel model = getModel ();
        Map<Integer, TableColumn> columnMap = 
            new HashMap<Integer, TableColumn>();
        for (Enumeration<TableColumn> en = columns.getColumns();
             en.hasMoreElements(); ) {
            TableColumn c = en.nextElement();
            columnMap.put(c.getModelIndex(), c);
        }

        for (int i = 0; i < model.getColumnCount(); ++i) {
            TableColumn c = columnMap.remove(i);
            if (c != null) {
                c.setHeaderValue(model.getColumnName(i));
            }
            else {
                c = new TableColumn (i);
                c.setPreferredWidth(100);
                addColumn (c);
            }
        }
		    
        // remove extra columns
        for (TableColumn c : columnMap.values()) {
            columns.removeColumn(c);
        }
        revalidate ();
        repaint ();
    }
}
