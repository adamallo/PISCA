/*
 * CNVSubstmodelParser.java
 *
 * By Diego Mallo
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package PISCA;

import dr.evolution.datatype.DataType;
import dr.evolution.datatype.GeneralDataType;
import PISCA.CNVSubstitutionModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.inference.model.Variable;
import dr.xml.*;

import java.util.logging.Logger;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @author Diego Mallo
 */
public class CNVSubstitutionModelParser extends AbstractXMLObjectParser {

	public static final String CNV_MODEL = "CNVModel";
	public static final String GAIN_RATE = "gain_rate";
    public static final String LOSS_RATE= "loss_rate";
    public static final String CONVERSION_RATE="conversion_rate";
    public static final String FREQUENCIES = "frequencies";

    public String getParserName() {
        return CNV_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Variable gainParam = (Variable) xo.getElementFirstChild(GAIN_RATE);
        Variable lossParam = (Variable) xo.getElementFirstChild(LOSS_RATE);
        Variable conversionParam = (Variable) xo.getElementFirstChild(CONVERSION_RATE);
        //FrequencyModel freqModel = (FrequencyModel) xo.getElementFirstChild(FrequencyModelParser.FREQUENCIES);
        XMLObject cxo = xo.getChild(FREQUENCIES);
        FrequencyModel freqModel = (FrequencyModel) cxo.getChild(FrequencyModel.class);
        DataType dataType = freqModel.getDataType();
        
        /*
        if (dataType != GeneralDataType.INSTANCE)
            throw new XMLParseException("Frequency model for the CNV substitution model must be coded as a generalDataType data type with 28 states: @ABCDEFGHIJKLMNOPQRSTUVWXYZ[ .");
         */
        Logger.getLogger("dr.evomodel").info("Creating CNV substitution model. Initial gain_rate = " +
                gainParam.getValue(0) + ", loss_rate = " + lossParam.getValue(0) + ", conversion_rate = " + conversionParam.getValue(0));

        return new CNVSubstitutionModel(dataType,gainParam,lossParam,conversionParam,freqModel);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element represents an instance of the CNV (Felsenstein, Kuhner, Maley, Mallo 2016) model of SGA evolution.";
    }

    public Class getReturnType() {
        return CNVSubstitutionModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
    		    new ElementRule(FREQUENCIES, FrequencyModel.class),
            //new ElementRule(FrequencyModelParser.FREQUENCIES,
             //       new XMLSyntaxRule[]{new ElementRule(FrequencyModel.class)}),
            new ElementRule(GAIN_RATE,
                    new XMLSyntaxRule[]{new ElementRule(Variable.class)}),
            new ElementRule(LOSS_RATE,
                    new XMLSyntaxRule[]{new ElementRule(Variable.class)}),
            new ElementRule(CONVERSION_RATE,
                    new XMLSyntaxRule[]{new ElementRule(Variable.class)})
    };
}
