import GFP_PV_PNN_Tools.Cell;
import GFP_PV_PNN_Tools.Tools;
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
import org.apache.commons.io.FilenameUtils;
import org.scijava.util.ArrayUtils;


/**
 * Detect DAPI nuclei, PV cells and PNN cells
 * Compute their colocalization
 * Detect DAPI and Gamma-H2AX foci in PNN/PV cells
 * @author Orion-CIRB
 */
public class GFP_PV_PNN implements PlugIn {
    
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
            // Write header in results file
            String header = "Image name\tCell label\tDAPI background\tNucleus vol (µm3)\tNucleus total int\t"
                    + "Nucleus total int corrected\tNucleus total Gamma-H2AX int\tNucleus total Gamma-H2AX int corrected\t" +
                    "PV background\tPV cell vol (µm3)\tPV cell total int\tPV cell total int corrected\t" +
                    "PNN background\tPNN cell vol (µm3)\tPNN cell total int\tPNN cell total int corrected\t" +
                    "Nb Gamma-H2AX foci\tGamma-H2AX foci total vol (µm3)\tGamma-H2AX foci total int\t" +
                    "Nb DAPI foci\tDAPI foci total vol (µm3)\tDAPI foci total int\n";
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
                tools.print("--- ANALYZING IMAGE " + rootName + " ------");
                reader.setId(f);
                
                ImporterOptions options = new ImporterOptions();
                options.setId(f);
                options.setSplitChannels(true);
                options.setQuiet(true);
                options.setColorMode(ImporterOptions.COLOR_MODE_GRAYSCALE);
                
                // Open DAPI channel
                tools.print("- Analyzing PV channel -");
                int indexCh = ArrayUtils.indexOf(chsName, channels[3]);
                ImagePlus imgDAPI = BF.openImagePlus(options)[indexCh];
                // Detect DAPI nuclei with CellPose
                System.out.println("Finding DAPI nuclei...");
                Objects3DIntPopulation dapiPop = tools.cellposeDetection(imgDAPI, true, tools.cellposeNucleiModel, 1, tools.cellposeNucleiDiameter, tools.cellposeNucleiCellThresh, 0.75, true, tools.minNucleusVol, tools.maxNucleusVol);
                System.out.println(dapiPop.getNbObjects() + " DAPI nuclei found");
                
                // Open PV channel
                tools.print("- Analyzing PV channel -");
                indexCh = ArrayUtils.indexOf(chsName, channels[0]);
                ImagePlus imgPV = BF.openImagePlus(options)[indexCh];
                // Detect PV cells with CellPose
                System.out.println("Finding PV cells....");
                Objects3DIntPopulation pvPop = tools.cellposeDetection(imgPV, true, tools.cellposePVModel, 1, tools.cellposePVDiameter, tools.cellposePVCellThresh, 0.25, false, tools.minCellVol, tools.maxCellVol);
                System.out.println(pvPop.getNbObjects() + " PV cells found");
                
                // Open PNN channel
                tools.print("- Analyzing PNN channel -");
                indexCh = ArrayUtils.indexOf(chsName, channels[1]);
                ImagePlus imgPNN = BF.openImagePlus(options)[indexCh];
                // Detect PNN cells with CellPose
                System.out.println("Finding PNN cells....");
                Objects3DIntPopulation pnnPop = tools.cellposeDetection(imgPNN, true, tools.cellposePNNModel, 1, tools.cellposePNNDiameter, tools.cellposePNNCellThresh, 0.25, false, tools.minCellVol, tools.maxCellVol);
                System.out.println(pnnPop.getNbObjects() + " PNN cells found");
                
                // Open GFP channel
                tools.print("- Analyzing GFP channel -");
                indexCh = ArrayUtils.indexOf(chsName, channels[2]);
                ImagePlus imgGFP = BF.openImagePlus(options)[indexCh];
                
                System.out.println("Colocalizing nuclei with PV and PNN cells....");
                ArrayList<Cell> cells = tools.colocalization(dapiPop, pvPop, pnnPop);
                System.out.println(cells.size() + " nuclei colocalized with a PV and/or a PNN cell");
                
                tools.print("- Measuring cells parameters -");
                tools.writeCellsParameters(cells, imgDAPI, imgPV, imgPNN, imgGFP);
               
                // Detect GFP foci in nuclei
                System.out.println("Finding GFP foci in each nucleus....");
                Objects3DIntPopulation gfpFociPop = tools.stardistFociInCellsPop(imgGFP, cells, "GFP", true);
                
                // Detect DAPI foci in nuclei
                System.out.println("Finding DAPI foci in each nucleus....");
                Objects3DIntPopulation dapiFociPop = tools.stardistFociInCellsPop(imgDAPI, cells, "DAPI", true);
                
                // Save image objects
                tools.print("- Saving results -");
                tools.drawResults(imgDAPI, imgPV, cells, gfpFociPop, dapiFociPop, rootName, outDirResults);
                
                // Write results
                for (Cell cell : cells) {
                    results.write(rootName+"\t"+cell.params.get("label")+"\t"+cell.params.get("dapiBg")+"\t"+cell.params.get("nucVol")+"\t"+cell.params.get("nucIntTot")+
                                  "\t"+cell.params.get("nucIntTotCorr")+"\t"+cell.params.get("nucGfpIntTot")+"\t"+cell.params.get("nucGfpIntTotCorr")+
                                  "\t"+cell.params.get("pvBg")+"\t"+cell.params.get("pvCellVol")+"\t"+cell.params.get("pvCellIntTot")+"\t"+cell.params.get("pvCellIntTotCorr")+
                                  "\t"+cell.params.get("pnnBg")+"\t"+cell.params.get("pnnCellVol")+"\t"+cell.params.get("pnnCellIntTot")+"\t"+cell.params.get("pnnCellIntTotCorr")+
                                  "\t"+cell.params.get("gfpFociNb")+"\t"+cell.params.get("gfpFociVolTot")+"\t"+cell.params.get("gfpFociIntTot")+
                                  "\t"+cell.params.get("dapiFociNb")+"\t"+cell.params.get("dapiFociVolTot")+"\t"+cell.params.get("dapiFociIntTot")+"\n");
                    results.flush();
                }
                
                tools.flush_close(imgDAPI);
                tools.flush_close(imgPNN);
                tools.flush_close(imgPV);
                tools.flush_close(imgGFP);
            }
        } catch (IOException | DependencyException | ServiceException | FormatException | io.scif.DependencyException  ex) {
            Logger.getLogger(GFP_PV_PNN.class.getName()).log(Level.SEVERE, null, ex);
        }
        tools.print("--- Process done ---");
    }    
}    
