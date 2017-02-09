/*
 * AscertainedCharactersSitePatternsParser.java
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

package pluginSCA;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.AscertainedSitePatterns;
import dr.evolution.alignment.PatternList;
import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.sequence.Sequence;
import dr.evolution.util.TaxonList;
import dr.util.Citable;
import dr.util.Identifiable;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.ContentRule;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

/**
 * AscertainedCharactersSitePatternsParser - based on AscertainedSitePatternsParser in order to add dummy patterns for given unsampled states
 *
 * @author Diego Mallo
 */
public class AscertainedCharactersSitePatternsParser extends AbstractXMLObjectParser  {
		public static final String ACPATTERNS="ascertainedCharacterPatterns";
		public static final String STATE="state";
		public static final String CODE="code";
	   
		public String getParserName() {
	        return ACPATTERNS;
	    }

	    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
	        SimpleAlignment alignment = (SimpleAlignment) xo.getChild(Alignment.class);
	        XMLObject xoc;
	        TaxonList taxa = null;
	        List<String> states = new ArrayList<String>();
	        
	        for (int i =0; i < xo.getChildCount(); i++) {
	            if (xo.getChild(i) instanceof XMLObject) {
	                xoc = (XMLObject)xo.getChild(i);
	                if (xoc.getName().equals(STATE)) {
	                    states.add(xoc.getStringAttribute(CODE));
	                } else {
	                    throw new XMLParseException("illegal element, " + xoc.getName() + ", in " + getParserName() + " element");
	                }
	            } else if (xo.getChild(i) instanceof Identifiable )  {
	            	    if(! (xo.getChild(i) instanceof Alignment))
	            	    {
	            	    		states.add(((Identifiable)xo.getChild(i)).getId());
	            	    }
	            } else
	            {
	                throw new XMLParseException("illegal element in " + getParserName() + " element");
	            }
	        }

	        int from = 0;
	        int to = alignment.getSiteCount()-1;


	        int startInclude = -1;//Unused
	        int stopInclude = -1;//Unused
	        int startExclude = alignment.getSiteCount();
	        int stopExclude = alignment.getSiteCount()+states.size();
	        
	        //Obtaining the dummy positions
	        StringBuilder outputStates=new StringBuilder(); //To print
	        StringBuilder seqToAppend=new StringBuilder(); //To append to the alignment
	        for(int i=0; i<states.size();++i) {
	        		String toappend = states.get(i);
	        		outputStates.append(toappend+",");
	        		seqToAppend.append(toappend);
	        }
	        outputStates.delete(outputStates.length()-1, outputStates.length());
	        
	        //New alignment with the dummy invariable patterns
	        SimpleAlignment newAlignment=new SimpleAlignment();
	        newAlignment.setDataType(alignment.getDataType());
	        for (int i = 0; i < alignment.getTaxonCount(); i++) {
	        		Sequence seq=alignment.getSequence(i);
	        		seq.appendSequenceString(seqToAppend.toString());
	        		newAlignment.addSequence(seq);
	        }

	        	//I should make sure that there is no invariable patterns with the given states
	        
	        Logger.getLogger("dr.evoxml").info("Creating ascertained site patterns '" + xo.getId() + "' including all positions of the alignment and correcting for the unsampled states " + outputStates.toString());


	        AscertainedSitePatterns patterns = new AscertainedSitePatterns(newAlignment, taxa,
	                from, to, 1,
	                startInclude, stopInclude,
	                startExclude, stopExclude);

	        Logger.getLogger("dr.evoxml").info("\tThere are " + patterns.getPatternCount() + " patterns in total, " + alignment.getSiteCount() + " original sites plus " + states.size() + " dummy sequences\n");
	        
	        Logger.getLogger("dr.evoxml").info("\tPlease cite:\n" + Citable.Utils.getCitationString(patterns));

	        return patterns;
	    }

	    public XMLSyntaxRule[] getSyntaxRules() {
	        return rules;
	    }

	    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
	            new ElementRule(Alignment.class),
	            new ElementRule(Identifiable.class, 0, Integer.MAX_VALUE),
	            new ContentRule("<state code=\"X\"/>")
	    };

	    public String getParserDescription() {
	        return "Generates a list of the unique site patterns (unique columns) in an alignment and the ascertainment correction of the unsampled states";
	    }
	    
	    public Class getReturnType() {
	        return PatternList.class;
	    }
}
