/*
 * BiallelicBinarySubstmodelParser.java
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
import PISCA.BiallelicBinarySubstitutionModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.xml.*;

import java.util.logging.Logger;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @author Diego Mallo
 */
public class BiallelicBinarySubstitutionModelParser extends AbstractXMLObjectParser {

	public static final String BiallelicBinary_MODEL = "BiallelicBinaryModel";
    public static final String DEMETHYLATION_RATE= "demethylation_rate";
    public static final String METHYLATION_RATE= "methylation_rate";
    public static final String FREQUENCIES = "frequencies";

    public String getParserName() {
        return BiallelicBinary_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

    	Variable demethylationParam=null;
    	Variable methylationParam=new Variable.D(1,1); //Not required. Default = 1.0
    	FrequencyModel freqModel=null;
    	DataType dataType=null;
    	
        for (int i = 0; i < xo.getChildCount(); i++) {
            XMLObject xoc = (XMLObject) xo.getChild(i);
            if (xoc.getName().equals(DEMETHYLATION_RATE)) {
                demethylationParam = (Variable) xoc.getChild(Parameter.class);
            } else if (xoc.getName().equals(METHYLATION_RATE)) {
                methylationParam = (Variable) xoc.getChild(Parameter.class);
            } else if (xoc.getName().equals(FREQUENCIES)) {
            	freqModel = (FrequencyModel) xoc.getChild(FrequencyModel.class);
            	dataType = freqModel.getDataType();
            }
            
        }

        Logger.getLogger("dr.evomodel").info("Creating BiallelicBinary substitution model. Initial methylation_rate = " + methylationParam.getValue(0) + "demethylation_rate = " +
        		demethylationParam.getValue(0));

        return new BiallelicBinarySubstitutionModel(dataType,methylationParam,demethylationParam,freqModel);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element represents an instance of the BiallelicBinary (Mallo, Malley 2022) model of  evolution.";
    }

    public Class getReturnType() {
        return BiallelicBinarySubstitutionModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
    		    new ElementRule(FREQUENCIES, FrequencyModel.class),
            //new ElementRule(FrequencyModelParser.FREQUENCIES,
             //       new XMLSyntaxRule[]{new ElementRule(FrequencyModel.class)}),
            new ElementRule(DEMETHYLATION_RATE,
                    new XMLSyntaxRule[]{new ElementRule(Variable.class)}),
            new ElementRule(METHYLATION_RATE,
                    new XMLSyntaxRule[]{new ElementRule(Variable.class)},true), //Not required
    };
}
