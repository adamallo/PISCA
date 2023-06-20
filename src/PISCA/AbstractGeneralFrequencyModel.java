package PISCA;

import dr.evolution.datatype.DataType;
import dr.evomodelxml.substmodel.FrequencyModelParser;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

/**
 * A model of equilibrium frequencies that does not need to have a frequencies parameter
 *
 * @author Diego M.
 */
public abstract class AbstractGeneralFrequencyModel extends AbstractModel {

    public AbstractGeneralFrequencyModel(DataType dataType) {

        super(FrequencyModelParser.FREQUENCY_MODEL);
        this.dataType = dataType;
    }

    public abstract void setFrequency(int i, double value);

    public abstract double getFrequency(int i);

    public abstract int getFrequencyCount();

    public double[] getFrequencies() {
        double[] frequencies = new double[getFrequencyCount()];
        for (int i = 0; i < frequencies.length; i++) {
            frequencies[i] = getFrequency(i);
        }
        return frequencies;
    }

    public double[] getCumulativeFrequencies() {
        double[] frequencies = getFrequencies();
        for (int i = 1; i < frequencies.length; i++) {
            frequencies[i] += frequencies[i - 1];
        }
        return frequencies;
    }

    public DataType getDataType() {
        return dataType;
    }

    // *****************************************************************
    // Interface Model
    // *****************************************************************

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        // no intermediates need recalculating....
    }

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        // no intermediates need recalculating....
    }

    protected void storeState() {
    } // no state apart from parameters to store

    protected void restoreState() {
    } // no state apart from parameters to restore

    protected void acceptState() {
    } // no state apart from parameters to accept

    public Element createElement(Document doc) {
        throw new RuntimeException("Not implemented!");
    }

    private DataType dataType = null;

}
