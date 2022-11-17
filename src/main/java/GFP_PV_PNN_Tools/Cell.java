package GFP_PV_PNN_Tools;

import java.util.HashMap;
import mcib3d.geom2.Object3DInt;

/**
 * @author Orion-CIRB
 */
public class Cell {
    
    // Nucleus
    private Object3DInt nucleus;
    // PV cell
    private Object3DInt pvCell;
    // PNN cell
    private Object3DInt pnnCell; 
    // Parameters
    public HashMap<String, Double> params;
    

    public Cell(Object3DInt nucleus) {
        this.nucleus = nucleus;
        this.pvCell = null;
        this.pnnCell = null;
        this.params = new HashMap<>();
    }
    
    
    public Object3DInt getNucleus() {
        return nucleus;
    }
    
    public Object3DInt getPvCell() {
        return pvCell;
    }
    
    public Object3DInt getPnnCell() {
        return pnnCell;
    }
    
    public void setLabel(double label) {
        params.put("label", label);
    }
    
    public void setNucParams(double dapiBg, double nucVol, double nucIntTot, double nucIntTotCorr, double gfpBg, double nucGfpIntTot, double nucGfpIntTotCorr) {
        params.put("dapiBg", dapiBg);
        params.put("nucVol", nucVol);
        params.put("nucIntTot", nucIntTot);
        params.put("nucIntTotCorr", nucIntTotCorr);
        params.put("gfpBg", gfpBg);
        params.put("nucGfpIntTot", nucGfpIntTot);
        params.put("nucGfpIntTotCorr", nucGfpIntTotCorr);
    }
    
    
    public void setPvCell(Object3DInt pvCell) {
        this.pvCell = pvCell;
    }
    
    public void setPvParams(double pvBg, Double pvCellVol, Double pvCellIntTot, Double pvCellIntTotCorr) {
        params.put("pvBg", pvBg);
        params.put("pvCellVol", pvCellVol);
        params.put("pvCellIntTot", pvCellIntTot);
        params.put("pvCellIntTotCorr", pvCellIntTotCorr);             
    }
    
    public void setPnnCell(Object3DInt pnnCell) {
        this.pnnCell = pnnCell;
    }
    
    public void setPnnParams(double pnnBg, Double pnnCellVol, Double pnnCellIntTot, Double pnnCellIntTotCorr) {
        params.put("pnnBg", pnnBg);
        params.put("pnnCellVol", pnnCellVol);
        params.put("pnnCellIntTot", pnnCellIntTot);
        params.put("pnnCellIntTotCorr", pnnCellIntTotCorr);   
    }
    
    public void setDapiFociParams(double dapiFociNb, double dapiFociVolTot, double dapiFociIntTot, double dapiFociIntTotCorr) {
        params.put("dapiFociNb", dapiFociNb);
        params.put("dapiFociVolTot", dapiFociVolTot);
        params.put("dapiFociIntTot", dapiFociIntTot);
        params.put("dapiFociIntTotCorr", dapiFociIntTotCorr);
    }
    
    public void setGfpFociParams(double gfpFociNb, double gfpFociVolTot, double gfpFociIntTot, double gfpFociIntTotCorr) {
        params.put("gfpFociNb", gfpFociNb);
        params.put("gfpFociVolTot", gfpFociVolTot);
        params.put("gfpFociIntTot", gfpFociIntTot);
        params.put("gfpFociIntTotCorr", gfpFociIntTotCorr);
    }
}
