/*
 * CNVDataType.java
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
 

package PISCA;

import dr.evolution.datatype.DataType;

import dr.util.Identifiable;

import java.util.*;

*//**
 * Implements a model for copy number variants
 *
 * @author Diego Mallo
 * @author Mary Kuhner
 * @author Jon Yamato
 * @author Joe Felsenstein 
 * @author Carlo Maley
 * @version $Id: CNVDataType.java,v 1.0 2016/02/10 20:25:56 $
 *//*
public class CNVDataType extends DataType implements Identifiable {

    public static final String CNV_DATA_TYPE = "cnvDataType";
    public static final String DESCRIPTION = CNV_DATA_TYPE;
    public static final String CNVM = "CNVM";
    public static final int CNV = 9;
    public static final int TYPE = CNV;
    public static final CNVDataType INSTANCE = new CNVDataType();

    public static final char [] CNV_CHARS={'A','G','M','S','Y','_','e','B','H','N','T','Z','`','f','C','I','O','U','[','a','g','D','J','P','V','\\','b','h','E','K','Q','W',']','c','i','F','L','R','X','^','d','j','G','M','S','Y','_','e','k',UNKNOWN_CHARACTER,GAP_CHARACTER};
*//**
 * Not developed further. Exploring the usage of the generalDatatype, since it would be difficult to add a new datatype from an external module. 
 *//*
    public static final int [] CNV_STATES={};
    
    protected CNVDataType() {
        for (int i = 0; i < stateCodes.length; i++) {
            State state = new State(i, stateCodes[i]);
            states.add(state);
            stateMap.put(stateCodes[i], state);
        }
        stateCount = states.size();

        this.ambiguousStateCount = 0;

    }

    *//**
     * Add an alias (a state code that represents a particular state).
     * Note that all this does is put an extra entry in the stateNumbers
     * array.
     *
     * @param alias a string that represents the state
     * @param code the state number
     *//*
    public void addAlias(String alias, String code) {
        State state =stateMap.get(code);
        if (state == null) {
            throw new IllegalArgumentException("DataType doesn't contain the state, " + code);
        }
        stateMap.put(alias, state);
    }

    *//**
     * Add an ambiguity (a state code that represents multiple states).
     *
     * @param code            a string that represents the state
     * @param ambiguousStates the set of states that this code refers to.
     *//*
    public void addAmbiguity(String code, String[] ambiguousStates) {

        int n = ambiguousStateCount + stateCount;

        int[] indices = new int[ambiguousStates.length];
        int i = 0;
        for (String stateCode : ambiguousStates) {
            State state =stateMap.get(stateCode);
            if (state == null) {
                throw new IllegalArgumentException("DataType doesn't contain the state, " + stateCode);
            }
            indices[i] = state.number;
            i++;
        }
        State state = new State(n, code, indices);
        states.add(state);
        ambiguousStateCount++;

        stateMap.put(code, state);
    }

    @Override
    public char[] getValidChars() {
        return null;
    }

    *//**
     * Get state corresponding to a code
     *
     * @param code string code
     * @return state
     *//*
    public int getState(String code) {
        if (code.equals("?")) {
            return getUnknownState();
        }
        if (!stateMap.containsKey(code)) {
            return -1;
        }
        return stateMap.get(code).number;
    }

    *//**
     * Override this function to cast to string codes...
     * @param c character
     *
     * @return
     *//*
    public int getState(char c) {
        return getState(String.valueOf(c));
    }
    *//**
     * Get state corresponding to an unknown
     *
     * @return state
     *//*
    public int getUnknownState() {
        return stateCount + ambiguousStateCount;
    }

    *//**
     * Get state corresponding to a gap
     *
     * @return state
     *//*
    public int getGapState() {
        return getUnknownState();
    }

    *//**
     * Get character corresponding to a given state
     *
     * @param state state
     * @return corresponding code
     *//*
    public String getCode(int state) {
        return states.get(state).code;
    }

    *//**
     * returns an array containing the non-ambiguous states
     * that this state represents.
     *//*
    public int[] getStates(int state) {

        return states.get(state).ambiguities;
    }

    *//**
     * returns an array containing the non-ambiguous states that this state represents.
     *//*
    public boolean[] getStateSet(int state) {

        if (state >= states.size()) {
            throw new IllegalArgumentException("invalid state index");
        }
        State s = states.get(state);

        boolean[] stateSet = new boolean[stateCount];
        for (int i = 0; i < stateCount; i++)
            stateSet[i] = false;

        for (int i = 0, n = s.ambiguities.length; i < n; i++) {
            stateSet[s.ambiguities[i]] = true;
        }

        return stateSet;
    }

    *//**
     * description of data type
     *
     * @return string describing the data type
     *//*
    public String getDescription() {
        if (id != null) {
            return id;
        } else {
            return DESCRIPTION;
        }
    }

    *//**
     * type of data type
     *
     * @return integer code for the data type
     *//*
    public int getType() {
        return TYPE;
    }

    // **************************************************************
    // Identifiable IMPLEMENTATION
    // **************************************************************

    private String id = null;

    public void setId(String id) {
        this.id = id;
    }

    public String getId() {
        return id;
    }

    private List<State> states = new ArrayList<State>();
    private Map<String, State> stateMap = new TreeMap<String, State>();

    private class State {
        int number;
        String code;

        int[] ambiguities;

        State(int number, String code) {
            this.number = number;
            this.code = code;
            this.ambiguities = new int[]{number};
        }

        State(int number, String code, int[] ambiguities) {
            this.number = number;
			this.code = code;
			this.ambiguities = ambiguities;
		}
	}

}
*/
