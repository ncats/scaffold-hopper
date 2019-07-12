package gov.nih.ncgc.rgroup;

import java.awt.Component;
import javax.swing.JTable;
import javax.swing.table.*;

public class ScaffoldDataCellRenderer extends DefaultTableCellRenderer {
    public ScaffoldDataCellRenderer () {
        setHorizontalAlignment (CENTER);
    }
    
    public Component getTableCellRendererComponent 
        (JTable table, Object value, boolean selected, 
         boolean focus, int row, int column) {
        if (value == null) {
            return super.getTableCellRendererComponent
                (table, value, selected, focus, row, column);
        }
        
        ScaffoldData data = (ScaffoldData)value;
        setText (String.format("<html>%1$.1f \u00B1%2$.2f",
                               data.getMean(), data.getStd()));
        setBackground (selected ? table.getSelectionBackground ()
                       : table.getBackground());
        
        return this;
    }
}
