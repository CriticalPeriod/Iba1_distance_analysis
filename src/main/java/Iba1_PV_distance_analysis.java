import ij.*;
import ij.plugin.PlugIn;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.logging.Level;
import java.util.logging.Logger;
import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.common.services.ServiceFactory;
import loci.formats.FormatException;
import loci.formats.meta.IMetadata;
import loci.formats.services.OMEXMLService;
import loci.plugins.BF;
import loci.plugins.util.ImageProcessorReader;
import loci.plugins.in.ImporterOptions;
import mcib3d.geom2.Objects3DIntPopulation;
import mcib3d.geom2.VoxelInt;
import mcib3d.geom2.measurements.MeasureCentroid;


import org.apache.commons.io.FilenameUtils;
import org.scijava.util.ArrayUtils;

import Iba1_PV_distance_analysis_Tools.Cell;
import Iba1_PV_distance_analysis_Tools.Tools;


/**
 * Detect DAPI nuclei, PV cells and PNN cells
 * Compute their colocalization
 * Detect DAPI and Gamma-H2AX foci in PNN/PV cells
 * @author Orion-CIRB
 */
public class Iba1_PV_distance_analysis implements PlugIn {
    
    Tools tools = new Tools();
    
    private String imageDir = "";
    public String outDirResults = "";
    public BufferedWriter results;
   
    
    public void run(String arg) {
        try {
            if ((!tools.checkInstalledModules()) || (!tools.checkStarDistModels())) {
                return;
            } 
            
            imageDir = IJ.getDirectory("Choose directory containing image files...");
            if (imageDir == null) {
                return;
            }   
            // Find images with extension
            String fileExt = tools.findImageType(new File(imageDir));
            ArrayList<String> imageFiles = tools.findImages(imageDir, fileExt);
            if (imageFiles == null) {
                IJ.showMessage("Error", "No images found with " + fileExt + " extension");
                return;
            }
            
            // Create output folder
            outDirResults = imageDir + File.separator+ "Results" + File.separator;
            File outDir = new File(outDirResults);
            if (!Files.exists(Paths.get(outDirResults))) {
                outDir.mkdir();
            }
            String header = "Image name\tCell label\tDAPI background\tNucleus vol (µm3)\tNucleus total int\tNucleus total int corrected\t" + 
            "Gamma-H2AX background\tNucleus total Gamma-H2AX int\tNucleus total Gamma-H2AX int corrected\t" +
            "PV background\tPV cell vol (µm3)\tPV cell total int\tPV cell total int corrected\t" +
            "Iba1 background\tIba1 cell vol (µm3)\tIba1 cell total int\tIba1 cell total int corrected\t" +
            "Mb cell vol (µm3)\tMb cell total int\tMb cell total int corrected\t" +
            "Nb dots \t dots total vol (µm3)\tdots foci total int\tdots foci total int corrected\t" +
            "Nb DAPI foci\tDAPI foci total vol (µm3)\tDAPI foci total int\tDAPI foci total int corrected\t" +
            "Centroid X\tCentroid Y\tCentroid Z\n";
            FileWriter fwResults = new FileWriter(outDirResults + "results.xls", false);
            results = new BufferedWriter(fwResults);
            results.write(header);
            results.flush();
            
            // Create OME-XML metadata store of the latest schema version
            ServiceFactory factory;
            factory = new ServiceFactory();
            OMEXMLService service = factory.getInstance(OMEXMLService.class);
            IMetadata meta = service.createOMEXMLMetadata();
            ImageProcessorReader reader = new ImageProcessorReader();
            reader.setMetadataStore(meta);
            reader.setId(imageFiles.get(0));
            
            // Find image calibration
            tools.cal = tools.findImageCalib(meta);
            
            // Find channel names
            String[] chsName = tools.findChannels(imageFiles.get(0), meta, reader);
            
            // Channels dialog
            String[] channels = tools.dialog(chsName);
            if (channels == null) {
                IJ.showStatus("Plugin canceled");
                return;
            }

            for (String f : imageFiles) {
                String rootName = FilenameUtils.getBaseName(f);
                tools.print("Hello Camille ! --- ANALYZING IMAGE " + rootName + " ------");
                reader.setId(f);
                
                ImporterOptions options = new ImporterOptions();
                options.setId(f);
                options.setSplitChannels(true);
                options.setQuiet(true);
                options.setColorMode(ImporterOptions.COLOR_MODE_GRAYSCALE);
                
                // Open DAPI channel
                tools.print("- Analyzing DAPI channel -");
                int indexCh = ArrayUtils.indexOf(chsName, channels[0]);
                ImagePlus imgDAPI = BF.openImagePlus(options)[indexCh];
                // Detect DAPI nuclei with CellPose
                System.out.println("Finding DAPI nuclei...");
                Objects3DIntPopulation dapiPop = tools.cellposeDetection(imgDAPI, true, tools.cellposeNucleiModel, 1, tools.cellposeNucleiDiameter, tools.cellposeNucleiCellThresh, tools.cellposeNucleiStitchThresh, true, tools.minNucleusVol, tools.maxNucleusVol, false,  0.4);
                System.out.println(dapiPop.getNbObjects() + " DAPI nuclei found");
                
                // Open PV channel
                tools.print("- Analyzing PV channel -");
                indexCh = ArrayUtils.indexOf(chsName, channels[3]);
                ImagePlus imgPV = BF.openImagePlus(options)[indexCh];
                // Detect PV cells with CellPose
                System.out.println("Finding PV cells....");
                Objects3DIntPopulation pvPop = tools.cellposeDetection(imgPV, true, tools.cellposePVModel, 1, tools.cellposePVCellDiameter, tools.cellposeCellThresh, tools.cellposeCellStitchThresh, false, tools.minCellVol, tools.maxCellVol, false, tools.PVflow_threshold);
                System.out.println(pvPop.getNbObjects() + " PV cells found");
                
                // Open Iba1 channel
                tools.print("- Analyzing Iba1 channel -");
                indexCh = ArrayUtils.indexOf(chsName, channels[2]);
                ImagePlus imgPNN = BF.openImagePlus(options)[indexCh];
                // Detect Iba1 cells with CellPose
                System.out.println("Finding Iba1 cells....");
                Objects3DIntPopulation pnnPop = tools.cellposeDetection(imgPNN, true, tools.cellposePNNModel, 1, tools.cellposeIba1CellDiameter, tools.cellposeCellThresh, tools.cellposeCellStitchThresh, false, tools.minNucleusVol, tools.maxCellVol, false, 0.6);
                System.out.println(pnnPop.getNbObjects() + " Iba1 cells found");
                
                // Open GFP channel
                tools.print("- Analyzing GFP channel -");
                indexCh = ArrayUtils.indexOf(chsName, channels[1]);
                ImagePlus imgGFP = BF.openImagePlus(options)[indexCh];
                
                System.out.println("Colocalizing nuclei with PV and PNN cells....");
                ArrayList<Cell> cells = tools.colocalization(dapiPop, pvPop, pnnPop, tools.nucleusPVColocThresh, tools.nucleusPNNColocThresh, tools.expansionType);
                System.out.println(cells.size() + " nuclei colocalized with a PV and/or a PNN cell");
                
                tools.print("- Measuring cells parameters -");
                tools.writeCellsParameters(cells, imgDAPI, imgPV, imgPNN, imgGFP);
                tools.flush_close(imgPNN);
                if (tools.detectFoci) {
                    // Detect GFP foci in membrane
                    System.out.println("Finding GFP foci in each membrane....");
                    Objects3DIntPopulation gfpFociPop = tools.stardistFociInMbPop(imgGFP, cells, "GFP", true);
                    tools.flush_close(imgGFP);
                    
                    // Detect DAPI foci in nuclei
                    System.out.println("Finding DAPI foci in each nucleus....");
                    Objects3DIntPopulation dapiFociPop = tools.stardistFociInCellsPop(imgDAPI, cells, "DAPI", true);
                    
                    // Save image objects with foci
                    tools.print("- Saving results with foci -");
                    tools.drawResults(imgDAPI, imgPV, cells, gfpFociPop, dapiFociPop, rootName, outDirResults);
                } else {
                    // Save image objects without foci
                    tools.print("- Saving results without foci -");
                    tools.drawResults_wofoci(imgDAPI, imgPV, cells, rootName, outDirResults);
                }
                
                // Write results
                for (Cell cell : cells) {
                    MeasureCentroid measureCentroid = new MeasureCentroid(cell.getNucleus());
                    VoxelInt centroid = measureCentroid.getCentroidRoundedAsVoxelInt();
                    cell.setCentroidCoordinates(centroid.getX(), centroid.getY(), centroid.getZ());
                    results.write(rootName + "\t" + cell.params.get("label") + "\t" + cell.params.get("dapiBg") + "\t" + cell.params.get("nucVol") + "\t" + cell.params.get("nucIntTot") +
                    "\t" + cell.params.get("nucIntTotCorr") + "\t" + cell.params.get("gfpBg") + "\t" + cell.params.get("nucGfpIntTot") + "\t" + cell.params.get("nucGfpIntTotCorr") +
                    "\t" + cell.params.get("pvBg") + "\t" + cell.params.get("pvCellVol") + "\t" + cell.params.get("pvCellIntTot") + "\t" + cell.params.get("pvCellIntTotCorr") +
                    "\t" + cell.params.get("pnnBg") + "\t" + cell.params.get("pnnCellVol") + "\t" + cell.params.get("pnnCellIntTot") + "\t" + cell.params.get("pnnCellIntTotCorr") +
                    "\t" + cell.params.get("membraneVol") + "\t" + cell.params.get("membraneIntTot") + "\t" + cell.params.get("membraneIntTotCorr") +
                    "\t" + cell.params.get("gfpFociNb") + "\t" + cell.params.get("gfpFociVolTot") + "\t" + cell.params.get("gfpFociIntTot") +
                    "\t" + cell.params.get("gfpFociIntTotCorr") + "\t" + cell.params.get("dapiFociNb") + "\t" + cell.params.get("dapiFociVolTot") + "\t" + cell.params.get("dapiFociIntTot") +
                    "\t" + cell.params.get("dapiFociIntTotCorr") + "\t" + cell.params.get("centroidX") + "\t" + cell.params.get("centroidY") + "\t" + cell.params.get("centroidZ") + "\n");
                results.flush();
                }
                tools.flush_close(imgDAPI);
                tools.flush_close(imgPV);
            }
            results.close();
        } catch (IOException | DependencyException | ServiceException | FormatException | io.scif.DependencyException  ex) {
            Logger.getLogger(Iba1_PV_distance_analysis.class.getName()).log(Level.SEVERE, null, ex);
        }
        tools.print("--- Process done ---");
    }    
    
}    
