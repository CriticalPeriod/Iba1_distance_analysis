package Iba1_PV_distance_analysis_Tools;


import Iba1_PV_distance_analysis_Tools.Cellpose.CellposeSegmentImgPlusAdvanced;
import Iba1_PV_distance_analysis_Tools.Cellpose.CellposeTaskSettings;
import Iba1_PV_distance_analysis_Tools.StardistOrion.StarDist2D;
import fiji.util.gui.GenericDialogPlus;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.Roi;
import ij.io.FileSaver;
import ij.measure.Calibration;
import ij.plugin.Duplicator;
import ij.plugin.RGBStackMerge;
import ij.plugin.ZProjector;
import ij.process.ImageProcessor;
import io.scif.DependencyException;
import java.awt.Color;
import java.awt.Font;
// import java.io.BufferedWriter;
import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Set;
import loci.common.services.ServiceException;
import loci.formats.FormatException;
import loci.formats.meta.IMetadata;
import loci.plugins.util.ImageProcessorReader;
import mcib3d.geom2.BoundingBox;
import mcib3d.geom2.Object3DInt;
import mcib3d.geom2.Object3DPlane;
import mcib3d.geom2.Objects3DIntPopulation;
import mcib3d.geom2.Objects3DIntPopulationComputation;
import mcib3d.geom2.VoxelInt;
import mcib3d.geom2.measurements.MeasureIntensity;
import mcib3d.geom2.measurements.MeasureVolume;
import mcib3d.geom2.measurementsPopulation.MeasurePopulationColocalisation;
import mcib3d.image3d.ImageHandler;
import net.haesleinhuepf.clij.clearcl.ClearCLBuffer;
import net.haesleinhuepf.clij2.CLIJ2;
import org.apache.commons.io.FilenameUtils;
import org.scijava.util.ArrayUtils;



/**
 * @author Orion-CIRB
 */
public class Tools {
    
    public String[] channelNames = {"DAPI","Gamma-H2AX","Iba1", "PV"};
    public Calibration cal;
    public double pixVol= 0;
    
    private CLIJ2 clij2 = CLIJ2.getInstance();
    
    // Nuclei and cells detection with Cellpose
    private String cellposeEnvDirPath = IJ.isWindows()? System.getProperty("user.home")+"\\miniconda3\\envs\\CellPose" : "/opt/anaconda3/envs/vglut/";
    public final String cellposeModelPath = IJ.isWindows()? System.getProperty("user.home")+"\\miniconda3\\envs\\CellPose\\models\\" : "/opt/anaconda3/envs/vglut/models/";
    private final boolean useGpu = true;
    public final String cellposeNucleiModel = "cyto2";
    public final int cellposeNucleiDiameter = 30;
    public int cellposeNucleiCellThresh = 0;
    public final double cellposeNucleiStitchThresh = 0.50;
    public double minNucleusVol = 40;
    public double maxNucleusVol = 2500;
    
    public final String cellposePVModel = cellposeModelPath+"cyto_PV1"; // need to add Cellpose models folder path if own model (for Windows only, not Linux)
    public final String cellposePNNModel = cellposeModelPath+"cyto2_Iba1_microglia";
    public double cellposeCellThresh = -3;
    public int cellposeIba1CellDiameter = 30;
    public int cellposePVCellDiameter = 40;
    public double PVflow_threshold = 0.6;
    public final double cellposeCellStitchThresh = 0.20;
    public double minCellVol = 40;
    public double maxCellVol = 5500;
    public double nucleusPVColocThresh = 0.2;
    public double nucleusPNNColocThresh = 0.15;
    private double stdIntTh = 0;
    public double expansionFactor = 0.4;
    public int expansionType = 0;
    
    // Foci detection with StarDist
    public boolean detectFoci = false; // Default value is false
    private final File modelsPath = new File(IJ.getDirectory("imagej")+File.separator+"models");
    private final String stardistFociModel = "fociRNA-1.2.zip";
    private Object syncObject = new Object();
    private final String stardistOutput = "Label Image";
    private final double stardistPercentileBottom = 0.2;
    private final double stardistPercentileTop = 99.8;
    private double stardistFociProbThresh = 0.3;
    private final double stardistFociOverlayThresh = 0.15;
    public double minFociVol = 0.05;
    public double maxFociVol = 15;
    
    
    /**
     * Display a message in the ImageJ console and status bar
     */
    public void print(String log) {
        System.out.println(log);
        IJ.showStatus(log);
    }
    
    
    /**
     * Check that needed modules are installed
     */
    public boolean checkInstalledModules() {
        ClassLoader loader = IJ.getClassLoader();
        try {
            loader.loadClass("net.haesleinhuepf.clij2.CLIJ2");
        } catch (ClassNotFoundException e) {
            IJ.log("CLIJ not installed, please install from update site");
            return false;
        }
        try {
            loader.loadClass("mcib3d.geom.Object3D");
        } catch (ClassNotFoundException e) {
            IJ.log("3D ImageJ Suite not installed, please install from update site");
            return false;
        }
        return true;
    }
    
    
    /**
     * Check that required StarDist models are present in Fiji models folder
     */
    public boolean checkStarDistModels() {
        FilenameFilter filter = (dir, name) -> name.endsWith(".zip");
        File[] modelList = modelsPath.listFiles(filter);
        int index = ArrayUtils.indexOf(modelList, new File(modelsPath+File.separator+stardistFociModel));
        if (index == -1) {
            IJ.showMessage("Error", stardistFociModel + " StarDist model not found, please add it in Fiji models folder");
            return false;
        }
        return true;
    }
    
    
     /**
     * Find image type
     */
    public String findImageType(File imagesFolder) {
        String ext = "";
        String[] files = imagesFolder.list();
        for (String name : files) {
            String fileExt = FilenameUtils.getExtension(name);
            switch (fileExt) {
                case "nd" :
                   ext = fileExt;
                   break;
                case "czi" :
                   ext = fileExt;
                   break;
                case "lif"  :
                    ext = fileExt;
                    break;
                case "ics" :
                    ext = fileExt;
                    break;
                case "ics2" :
                    ext = fileExt;
                    break;
                case "tif" :
                    ext = fileExt;
                    break;
            }
        }
        return(ext);
    }
    
    
    /**
     * Find images in folder
     */
    public ArrayList<String> findImages(String imagesFolder, String imageExt) {
        File inDir = new File(imagesFolder);
        String[] files = inDir.list();
        if (files == null) {
            System.out.println("No image found in " + imagesFolder);
            return null;
        }
        ArrayList<String> images = new ArrayList();
        for (String f : files) {
            // Find images with extension
            String fileExt = FilenameUtils.getExtension(f);
            if (fileExt.equals(imageExt) && !f.startsWith("."))
                images.add(imagesFolder + File.separator + f);
        }
        Collections.sort(images);
        return(images);
    }
    
    
     /**
     * Find image calibration
     */
    public Calibration findImageCalib(IMetadata meta) {
        cal = new Calibration();
        cal.pixelWidth = meta.getPixelsPhysicalSizeX(0).value().doubleValue();
        cal.pixelHeight = cal.pixelWidth;
        if (meta.getPixelsPhysicalSizeZ(0) != null)
            cal.pixelDepth = meta.getPixelsPhysicalSizeZ(0).value().doubleValue();
        else
            cal.pixelDepth = 1;
        cal.setUnit("microns");
        System.out.println("XY calibration = " + cal.pixelWidth + ", Z calibration = " + cal.pixelDepth);
        return(cal);
    }
     /**
     * Find channels name
     * @param imageName
     * @return 
     * @throws loci.common.services.DependencyException
     * @throws loci.common.services.ServiceException
     * @throws loci.formats.FormatException
     * @throws java.io.IOException
     */
    public String[] findChannels (String imageName, IMetadata meta, ImageProcessorReader reader) throws DependencyException, ServiceException, FormatException, IOException {
        ArrayList<String> channels = new ArrayList<>();
        int chs = reader.getSizeC();
        String imageExt =  FilenameUtils.getExtension(imageName);
        switch (imageExt) {
            case "nd" :
                for (int n = 0; n < chs; n++) 
                {
                    if (meta.getChannelID(0, n) == null)
                        channels.add(Integer.toString(n));
                    else 
                        channels.add(meta.getChannelName(0, n).toString());
                }
                break;
            case "nd2" :
                for (int n = 0; n < chs; n++) 
                {
                    if (meta.getChannelID(0, n) == null)
                        channels.add(Integer.toString(n));
                    else 
                        channels.add(meta.getChannelName(0, n).toString());
                }
                break;
            case "lif" :
                for (int n = 0; n < chs; n++) 
                    if (meta.getChannelID(0, n) == null || meta.getChannelName(0, n) == null)
                        channels.add(Integer.toString(n));
                    else 
                        channels.add(meta.getChannelName(0, n).toString());
                break;
            case "czi" :
                for (int n = 0; n < chs; n++) 
                    if (meta.getChannelID(0, n) == null)
                        channels.add(Integer.toString(n));
                    else 
                        channels.add(meta.getChannelFluor(0, n).toString());
                break;
            case "ics" :
                for (int n = 0; n < chs; n++) 
                    if (meta.getChannelID(0, n) == null)
                        channels.add(Integer.toString(n));
                    else 
                        channels.add(meta.getChannelExcitationWavelength(0, n).value().toString());
                break; 
            case "ics2" :
                for (int n = 0; n < chs; n++) 
                    if (meta.getChannelID(0, n) == null)
                        channels.add(Integer.toString(n));
                    else 
                        channels.add(meta.getChannelExcitationWavelength(0, n).value().toString());
                break;        
            default :
                for (int n = 0; n < chs; n++)
                    channels.add(Integer.toString(n));

        }
        channels.add("None");
        return(channels.toArray(new String[channels.size()]));         
    }
    
    /**
     * Generate dialog box
     */
    public String[] dialog(String[] chs) {
        GenericDialogPlus gd = new GenericDialogPlus("Parameters");
        
        gd.addMessage("Channels", Font.getFont("Monospace"), Color.blue);
        int index = 0;
        for (String chNames: channelNames) {
            gd.addChoice(chNames + ": ", chs, chs[index]);
            index++;
        }
        
        gd.addMessage("Nuclei detection", Font.getFont("Monospace"), Color.blue);
        gd.addDirectoryField("Cellpose environment path: ", cellposeEnvDirPath);
        gd.addNumericField("Min nucleus volume (µm3): ", minNucleusVol);
        gd.addNumericField("Max nucleus volume (µm3): ", maxNucleusVol);
    
        gd.addMessage("Cells detection", Font.getFont("Monospace"), Color.blue);
        gd.addNumericField("Min cell volume (µm3): ", minCellVol);
        gd.addNumericField("Max cell volume (µm3): ", maxCellVol);
        gd.addNumericField("PNN intensity STD threshold: ", stdIntTh);
        gd.addNumericField("PV prob threshold: (-6 to 6) ", cellposeCellThresh);
        gd.addNumericField("PV intensity flow: ) ", PVflow_threshold);
        gd.addNumericField("Expansion Factor : ) ", expansionFactor);
        gd.addMessage("Expansion Type", Font.getFont("Monospace"), Color.blue);
        String[] expansionTypes = {"Nucleus", "PV", "PNN"};
        gd.addChoice("Expansion Type: ", expansionTypes, expansionTypes[0]);
        
        
        gd.addMessage("Foci detection", Font.getFont("Monospace"), Color.blue);
        gd.addCheckbox("Detect foci: ", detectFoci);
        gd.addNumericField("StarDist probability threshold", stardistFociProbThresh);
        gd.addNumericField("Min foci volume (µm3): ", minFociVol);
        gd.addNumericField("Max foci volume (µm3): ", maxFociVol);
        
        gd.addMessage("Image calibration", Font.getFont("Monospace"), Color.blue);
        gd.addNumericField("XY pixel size (µm): ", cal.pixelWidth);
        gd.addNumericField("Z pixel depth (µm):", cal.pixelDepth);
        gd.showDialog();
        
        String[] chChoices = new String[channelNames.length];
        for (int n = 0; n < chChoices.length; n++) 
            chChoices[n] = gd.getNextChoice();
        if(gd.wasCanceled())
            chChoices = null;
        
        cellposeEnvDirPath = gd.getNextString();
        minNucleusVol = gd.getNextNumber();
        maxNucleusVol = gd.getNextNumber();
        minCellVol = gd.getNextNumber();
        maxCellVol = gd.getNextNumber();
        stdIntTh = gd.getNextNumber();
        cellposeCellThresh = gd.getNextNumber();
        PVflow_threshold = gd.getNextNumber();
        expansionFactor = gd.getNextNumber();
        expansionType = gd.getNextChoiceIndex();
        detectFoci = gd.getNextBoolean();
        stardistFociProbThresh = gd.getNextNumber();
        minFociVol = gd.getNextNumber();
        maxFociVol = gd.getNextNumber();
        cal.pixelWidth = cal.pixelHeight = gd.getNextNumber();
        cal.pixelDepth = gd.getNextNumber();
        pixVol = cal.pixelWidth*cal.pixelWidth*cal.pixelDepth;
        
        return(chChoices);
    }
    
    
    /**
     * Flush and close an image
     */
    public void flush_close(ImagePlus img) {
        img.flush();
        img.close();
    }
    
    
    /**
     * Look for all 3D cells in a Z-stack: 
     * - apply CellPose in 2D slice by slice 
     * - let CellPose reconstruct cells in 3D using the stitch threshold parameters
     */
    public Objects3DIntPopulation cellposeDetection(ImagePlus img, boolean resize, String cellposeModel, int channel, int diameter, double cellProbThresh, double stitchThreshold, boolean removeOutliers, double volMin, double volMax, boolean intFilter, double PVflow_threshold) throws IOException{
        ImagePlus imgResized;
        if (resize) {
            float resizeFactor = 0.25f;
            imgResized = img.resize((int)(img.getWidth()*resizeFactor), (int)(img.getHeight()*resizeFactor), 1, "none");
        } else {
            imgResized = new Duplicator().run(img);
        }
        if (removeOutliers)
            IJ.run(imgResized, "Remove Outliers", "block_radius_x=40 block_radius_y=40 standard_deviations=1 stack");

        // Define CellPose settings
        CellposeTaskSettings settings = new CellposeTaskSettings(cellposeModel, channel, diameter, cellposeEnvDirPath);
        settings.setCellProbTh(cellProbThresh);
        settings.setStitchThreshold(stitchThreshold);
        settings.useGpu(useGpu);
        settings.setFlowTh(PVflow_threshold);
       
        // Run CellPose
        CellposeSegmentImgPlusAdvanced cellpose = new CellposeSegmentImgPlusAdvanced(settings, imgResized);
        ImagePlus imgOut = cellpose.run();
        if(resize) imgOut = imgOut.resize(img.getWidth(), img.getHeight(), "none");
        imgOut.setCalibration(cal);
       
        // Get cells as a population of objects and filter them
        ImageHandler imgH = ImageHandler.wrap(imgOut);
        Objects3DIntPopulation pop = new Objects3DIntPopulation(imgH);
        Objects3DIntPopulation popExcludeBorders = new Objects3DIntPopulationComputation(pop).getExcludeBorders(ImageHandler.wrap(img), false);
        Objects3DIntPopulationComputation popComputation = new Objects3DIntPopulationComputation(popExcludeBorders);
        Objects3DIntPopulation popFilter = popComputation.getFilterSize(volMin/pixVol, volMax/pixVol);
        Objects3DIntPopulation popFilterZ = zFilterPop(popFilter);
        if(intFilter)
            filterCellsByIntensitySTD(popFilterZ, img, stdIntTh);
        popFilterZ.resetLabels();
        System.out.println(popFilterZ.getNbObjects() + " detections remaining after filtering (" + (pop.getNbObjects()-popFilterZ.getNbObjects()) + " filtered out)");
               
        flush_close(imgOut);
        imgH.closeImagePlus();
        return(popFilterZ);
    }
    
    /**
     * Remove object with size < min and size > max in microns
     * @param pop
     * @param min
     * @param max
     */
    public void popFilterSize(Objects3DIntPopulation pop, double min, double max) {
        pop.getObjects3DInt().removeIf(p -> (new MeasureVolume(p).getVolumeUnit() < min) || (new MeasureVolume(p).getVolumeUnit() > max));
        pop.resetLabels();
    }
    
    /*
     * Remove objects present in only one z slice from population 
     */
    public Objects3DIntPopulation zFilterPop (Objects3DIntPopulation pop) {
        Objects3DIntPopulation popZ = new Objects3DIntPopulation();
        for (Object3DInt obj : pop.getObjects3DInt()) {
            if (obj.getBoundingBox().zmax != obj.getBoundingBox().zmin)
                popZ.addObject(obj);
        }
        return popZ;
    }
    
   
    /**
     * Remove cells if intensity STD is less than a certain threshold
     */
    public void filterCellsByIntensitySTD(Objects3DIntPopulation pop, ImagePlus img, double intStdTh) {
        ImageHandler imh = ImageHandler.wrap(img);
        pop.getObjects3DInt().removeIf(p -> (new MeasureIntensity(p, imh).getValueMeasurement(MeasureIntensity.INTENSITY_SD) <= intStdTh)); 
    }
        
    
    /**
     * Find cells colocalizing with a nucleus
     */
    public ArrayList<Cell> colocalization(Objects3DIntPopulation nucleiPop, Objects3DIntPopulation cellPop1, Objects3DIntPopulation cellPop2, double colocThresh1, double colocThresh2, int expansionType) {
        ArrayList<Cell> cells = new ArrayList<Cell>();
        if (nucleiPop.getNbObjects() > 0) {
            MeasurePopulationColocalisation coloc1 = new MeasurePopulationColocalisation(nucleiPop, cellPop1);
            MeasurePopulationColocalisation coloc2 = new MeasurePopulationColocalisation(nucleiPop, cellPop2);
            float label = 1;
            
            for (Object3DInt nucleus: nucleiPop.getObjects3DInt()) {
                Cell cell = new Cell(nucleus);
                boolean coloc = false;
             
                for (Object3DInt c1: cellPop1.getObjects3DInt()) {
                    double colocVal = coloc1.getValueObjectsPair(nucleus, c1);
                    if (colocVal > colocThresh1*nucleus.size()) {
                        cell.setPvCell(c1);
                        cellPop1.removeObject(c1);
                        coloc = true;
                        break;
                    }
                }
                
                for (Object3DInt c2: cellPop2.getObjects3DInt()) {
                    double colocVal = coloc2.getValueObjectsPair(nucleus, c2);
                    if (colocVal > colocThresh2*nucleus.size()) {
                        cell.setPnnCell(c2);
                        cellPop2.removeObject(c2);
                        coloc = true;
                        break;
                    }
                }
                
                if (coloc) {
                    // Extend the specified object by 40% to create membrane
                    Object3DInt expansionObject = null;
                    switch (expansionType) {
                        case 0: // Nucleus
                            expansionObject = nucleus;
                            break;
                        case 1: // PV
                            expansionObject = cell.getPvCell();
                            break;
                        case 2: // PNN
                            expansionObject = cell.getPnnCell();
                            break;
                    }
    
                    if (expansionObject != null) {
                        Object3DInt membrane = extendObject_c(expansionObject, expansionFactor);
                        cell.setMembrane(membrane);

                    }
    
                    cell.setLabel(label);
                    cells.add(cell);
                    label++;
                }

            }
        }
        return(cells);
    }
    private Object3DInt extendObject_c(Object3DInt object, double expansionFactor) {
        // Create a new Object3DInt for the expanded object (the shell)
        Object3DInt shellObject = new Object3DInt();
        
        // Get the bounding box of the object
        BoundingBox boundingBox = object.getBoundingBox();
        int xMin = boundingBox.xmin;
        int xMax = boundingBox.xmax;
        int yMin = boundingBox.ymin;
        int yMax = boundingBox.ymax;
        int zMin = boundingBox.zmin;
        int zMax = boundingBox.zmax;
        
        // Calculate the expansion distance in voxels (using the max dimension as reference)
        int xSize = xMax - xMin;
        int ySize = yMax - yMin;
        int zSize = zMax - zMin;
        int maxSize = Math.max(Math.max(xSize, ySize), zSize);
        int expansionDistance = (int)Math.ceil(maxSize * expansionFactor / 2);
        
        // Extend the bounding box for the dilation
        xMin -= expansionDistance;
        xMax += expansionDistance;
        yMin -= expansionDistance;
        yMax += expansionDistance;
        zMin -= expansionDistance;
        zMax += expansionDistance;
        
        // Create a set to store the original object voxels for quick lookup
        Set<String> originalVoxelKeys = new HashSet<>();
        for (Object3DPlane plane : object.getObject3DPlanes()) {
            for (VoxelInt voxel : plane.getVoxels()) {
                originalVoxelKeys.add(voxel.getX() + "_" + voxel.getY() + "_" + voxel.getZ());
            }
        }
        
        // Calculate the center of the object
        double[] centroid = calculateCentroid(object);
        
        // For each voxel in the extended bounding box, check if it should be part of the shell
        for (int z = zMin; z <= zMax; z++) {
            for (int y = yMin; y <= yMax; y++) {
                for (int x = xMin; x <= xMax; x++) {
                    // Calculate distance from centroid to current point
                    double dx = x - centroid[0];
                    double dy = y - centroid[1];
                    double dz = z - centroid[2];
                    // double distance = Math.sqrt(dx*dx + dy*dy + dz*dz);
                    
                    // Calculate distance from centroid to current point in original object scale
                    // double originalDistance = distance / (1 + expansionFactor);
                    
                    // Calculate the coordinates this would map to in the original object
                    int origX = (int)Math.round(centroid[0] + dx * (1 / (1 + expansionFactor)));
                    int origY = (int)Math.round(centroid[1] + dy * (1 / (1 + expansionFactor)));
                    int origZ = (int)Math.round(centroid[2] + dz * (1 / (1 + expansionFactor)));
                    
                    String newKey = x + "_" + y + "_" + z;
                    String origKey = origX + "_" + origY + "_" + origZ;
                    
                    // If this point maps to the original object but is not in the original object,
                    // then it's part of the shell
                    if (originalVoxelKeys.contains(origKey) && !originalVoxelKeys.contains(newKey)) {
                        shellObject.addVoxel(new VoxelInt(x, y, z, object.getLabel()));
                    }
                }
            }
        }
        shellObject.setVoxelSizeXY(cal.pixelWidth);
        shellObject.setVoxelSizeZ(cal.pixelDepth);
        
        return shellObject;
    }
    // Old version of the extension. Contains a mistake, fixed in the above version. 
    // private Object3DInt extendObject(Object3DInt object, double expansionFactor) {
    //     // Calculate the centroid of the object
    //     double[] centroid = calculateCentroid(object);
    
    //     // Create a new Object3DInt for the expanded object
    //     Object3DInt expandedObject = new Object3DInt();
    
    //     // Get the bounding box of the object to determine the z-range
    //     BoundingBox boundingBox = object.getBoundingBox();
    //     int zMin = boundingBox.zmin;
    //     int zMax = boundingBox.zmax;
    
    //     // Create a set to store the original object voxels for quick lookup
    //     Set<VoxelInt> originalVoxels = new HashSet<>();
    //     for (Object3DPlane plane : object.getObject3DPlanes()) {
    //         originalVoxels.addAll(plane.getVoxels());
    //     }
    
    //     // Iterate through the planes and their voxels
    //     for (Object3DPlane plane : object.getObject3DPlanes()) {
    //         for (VoxelInt voxel : plane.getVoxels()) {
    //             int x = voxel.getX();
    //             int y = voxel.getY();
    //             int z = voxel.getZ();
    
    //             // Calculate the new coordinates based on the expansion factor
    //             int newX = (int) Math.round(centroid[0] + (x - centroid[0]) * (1 + expansionFactor));
    //             int newY = (int) Math.round(centroid[1] + (y - centroid[1]) * (1 + expansionFactor));
    //             int newZ = (int) Math.round(centroid[2] + (z - centroid[2]) * (1 + expansionFactor));
    
    //             // Ensure the new z-coordinate is within the original z-range
    //             newZ = Math.max(zMin, Math.min(zMax, newZ));
    
    //             // Create the new voxel
    //             VoxelInt newVoxel = new VoxelInt(newX, newY, newZ, voxel.getValue());
    
    //             // Add the expanded voxel to the new object if it's not part of the original object
    //             if (!originalVoxels.contains(newVoxel)) {
    //                 expandedObject.addVoxel(newVoxel);
    //             }
    //         }
    //     }
    
    //     return expandedObject;
    // }
        
    private double[] calculateCentroid(Object3DInt object) {
        double sumX = 0;
        double sumY = 0;
        double sumZ = 0;
        int count = 0;
    
        for (Object3DPlane plane : object.getObject3DPlanes()) {
            for (VoxelInt voxel : plane.getVoxels()) {
                sumX += voxel.getX();
                sumY += voxel.getY();
                sumZ += voxel.getZ();
                count++;
            }
        }
    
        return new double[]{sumX / count, sumY / count, sumZ / count};
    }

    /**
     * Compute and save PV and PNN cells parameters
     */
    public void writeCellsParameters(ArrayList<Cell> cells, ImagePlus imgDAPI, ImagePlus imgPV, ImagePlus imgPNN, ImagePlus imgGFP) {
        double dapiBg = findBackground(imgDAPI, "DAPI");
        double pvBg = findBackground(imgPV, "PV");
        double pnnBg = findBackground(imgPNN, "PNN");
        double gfpBg = findBackground(imgGFP, "GFP");
        
        for (Cell cell : cells) {
            float label = cell.params.get("label").floatValue();
            
            // Nucleus
            Object3DInt nucleus = cell.getNucleus();
            nucleus.setLabel(label);
            double nucVol = new MeasureVolume(nucleus).getVolumeUnit();
            double nucIntTot = new MeasureIntensity(nucleus, ImageHandler.wrap(imgDAPI)).getValueMeasurement(MeasureIntensity.INTENSITY_SUM);
            double nucGFPIntTot = new MeasureIntensity(nucleus, ImageHandler.wrap(imgGFP)).getValueMeasurement(MeasureIntensity.INTENSITY_SUM);
            cell.setNucParams(dapiBg, nucVol, nucIntTot, nucIntTot-dapiBg*nucVol, gfpBg, nucGFPIntTot, nucGFPIntTot-gfpBg*nucVol);

            // PV cell
            Object3DInt pvCell = cell.getPvCell();
            if(pvCell != null) {
                pvCell.setLabel(label);
                double pvCellVol = new MeasureVolume(pvCell).getVolumeUnit();
                double pvCellIntTot = new MeasureIntensity(pvCell, ImageHandler.wrap(imgPV)).getValueMeasurement(MeasureIntensity.INTENSITY_SUM);
                cell.setPvParams(pvBg, pvCellVol, pvCellIntTot, pvCellIntTot-pvBg*pvCellVol);
            } else {
                cell.setPvParams(pvBg, null, null, null);
            }
            
            // PNN cell
            Object3DInt pnnCell = cell.getPnnCell();
            if(pnnCell != null) {
                pnnCell.setLabel(label);
                double pnnCellVol = new MeasureVolume(pnnCell).getVolumeUnit();
                double pnnCellIntTot = new MeasureIntensity(pnnCell, ImageHandler.wrap(imgPNN)).getValueMeasurement(MeasureIntensity.INTENSITY_SUM);
                cell.setPnnParams(pnnBg, pnnCellVol, pnnCellIntTot, pnnCellIntTot-pnnBg*pnnCellVol);
            } else {
                cell.setPnnParams(pnnBg, null, null, null);
            }
            // Membrane
            Object3DInt membrane = cell.getMembrane();
            if(membrane != null) {
                membrane.setLabel(label);
                double membraneVol = new MeasureVolume(membrane).getVolumeUnit();
                double membraneIntTot = new MeasureIntensity(membrane, ImageHandler.wrap(imgGFP)).getValueMeasurement(MeasureIntensity.INTENSITY_SUM);
                cell.setMembraneParams(membraneVol, membraneIntTot, membraneIntTot-gfpBg*membraneVol);
            }
        }
    }
    
    
    /**
     * Do Z projection
     * @param img
     * @param projection parameter
     */
    public ImagePlus doZProjection(ImagePlus img, int param) {
        ZProjector zproject = new ZProjector();
        zproject.setMethod(param);
        zproject.setStartSlice(1);
        zproject.setStopSlice(img.getNSlices());
        zproject.setImage(img);
        zproject.doProjection();
       return(zproject.getProjection());
    }

    /**
     * Find background image intensity:
     * Z projection over min intensity + read median intensity
     * @param img
     */
    public double findBackground(ImagePlus img, String channelName) {
      ImagePlus imgProj = doZProjection(img, ZProjector.MIN_METHOD);
      ImageProcessor imp = imgProj.getProcessor();
      double bg = imp.getStatistics().median;
      System.out.println(channelName + " background (median of the min projection) = " + bg);
      flush_close(imgProj);
      return(bg);
    }
    
    
    /** 
    * For each nucleus find foci
    * return foci pop cell population
    */
    public Objects3DIntPopulation stardistFociInCellsPop(ImagePlus img, ArrayList<Cell> cells, String fociType, boolean resize) throws IOException{
        float fociIndex = 1;
        double resizeFactor = 0.5;
        Objects3DIntPopulation allFociPop = new Objects3DIntPopulation();
        for (Cell cell: cells) {
            Object3DInt nuc = cell.getNucleus();
            
            // Crop image around nucleus
            BoundingBox box = nuc.getBoundingBox();
            Roi roiBox = new Roi(box.xmin, box.ymin, box.xmax-box.xmin, box.ymax-box.ymin);
            img.setRoi(roiBox);
            img.updateAndDraw();
            ImagePlus imgNuc = new Duplicator().run(img, box.zmin+1, box.zmax+1);
            imgNuc.deleteRoi();
            imgNuc.updateAndDraw();
            
            // Median filter, downscaling and gaussian filter
            ImagePlus imgM = median_filter(imgNuc, 1, 1);
            ImagePlus imgS = (resize) ? imgM.resize((int)(resizeFactor*imgNuc.getWidth()), (int)(resizeFactor*imgNuc.getHeight()), 1, "average") : imgM.duplicate();
            flush_close(imgM);
            ImagePlus imgG = gaussian_filter(imgS,1, 1);
            flush_close(imgS);

            // StarDist
            File starDistModelFile = new File(modelsPath+File.separator+stardistFociModel);
            StarDist2D star = new StarDist2D(syncObject, starDistModelFile);
            star.loadInput(imgG);
            star.setParams(stardistPercentileBottom, stardistPercentileTop, stardistFociProbThresh, stardistFociOverlayThresh, stardistOutput);
            star.run();
            flush_close(imgG);

            // Label foci in 3D
            ImagePlus imgLabels = star.associateLabels();
            if (resize) imgLabels = imgLabels.resize(imgNuc.getWidth(), imgNuc.getHeight(), 1, "none");
            
            flush_close(imgNuc);
            imgLabels.setCalibration(cal);
            Objects3DIntPopulation fociPop = new Objects3DIntPopulation(ImageHandler.wrap(imgLabels));
            popFilterSize(fociPop, minFociVol, maxFociVol);
            flush_close(imgLabels);
            
            // Find foci in nucleus
            fociPop.translateObjects(box.xmin, box.ymin, box.zmin);
            Objects3DIntPopulation fociColocPop = findFociInCell(nuc, fociPop);
            System.out.println(fociColocPop.getNbObjects() + " " + fociType + " foci found in nucleus " + nuc.getLabel());
            
            for (Object3DInt foci: fociColocPop.getObjects3DInt()) {
                foci.setLabel(fociIndex);
                fociIndex++;
                foci.setType(cell.params.get("label").intValue());
                allFociPop.addObject(foci);
            }
            
            writeFociParameters(cell, fociColocPop, img, fociType);
        }
        return(allFociPop);
    }
    
    
    /**
     * Median filter using CLIJ2
     */ 
    public ImagePlus median_filter(ImagePlus img, double sizeXY, double sizeZ) {
       ClearCLBuffer imgCL = clij2.push(img);
       ClearCLBuffer imgCLMed = clij2.create(imgCL);
       clij2.median3DBox(imgCL, imgCLMed, sizeXY, sizeXY, sizeZ);
       clij2.release(imgCL);
       ImagePlus imgMed = clij2.pull(imgCLMed);
       clij2.release(imgCLMed);
       return(imgMed);
    }
    
    /**
     * Guassian filter using CLIJ2
     * @param img
     * @param sizeXY
     * @param sizeZ
     * @return 
     */ 
    public ImagePlus gaussian_filter(ImagePlus img, double sizeXY, double sizeZ) {
       ClearCLBuffer imgCL = clij2.push(img); 
       ClearCLBuffer imgCLGauss = clij2.create(imgCL);
       clij2.gaussianBlur3D(imgCL, imgCLGauss, sizeXY, sizeXY, sizeZ);
       clij2.release(imgCL);
       ImagePlus imgGauss = clij2.pull(imgCLGauss);
        clij2.release(imgCLGauss);
       return(imgGauss);
    } 
    
    
    /**
     * Find dots population colocalizing with a cell objet
     */
    public Objects3DIntPopulation findFociInCell(Object3DInt cellObj, Objects3DIntPopulation dotsPop) {
        Objects3DIntPopulation cellPop = new Objects3DIntPopulation();
        cellPop.addObject(cellObj);
        Objects3DIntPopulation colocPop = new Objects3DIntPopulation();
        if (dotsPop.getNbObjects() > 0) {
            MeasurePopulationColocalisation coloc = new MeasurePopulationColocalisation(cellPop, dotsPop);
            for (Object3DInt dot: dotsPop.getObjects3DInt()) {
                    double colocVal = coloc.getValueObjectsPair(cellObj, dot);
                    if (colocVal > 0.5*dot.size()) {
                        colocPop.addObject(dot);
                    }
            }
        }
        return(colocPop);
    }

    
    /**
     * Compute foci parameters and save them in corresponding Nucleus
     */
    public void writeFociParameters(Cell cell, Objects3DIntPopulation fociColocPop, ImagePlus img, String fociType) {
        int fociNb = fociColocPop.getNbObjects();
        double fociVol = 0;
        double fociInt = 0;
        ImageHandler imh = ImageHandler.wrap(img.duplicate());
        for (Object3DInt obj: fociColocPop.getObjects3DInt()) {
            fociVol += new MeasureVolume(obj).getVolumeUnit();
            fociInt += new MeasureIntensity(obj, imh).getValueMeasurement(MeasureIntensity.INTENSITY_SUM);
        }
        
        switch (fociType) {
            case "DAPI" :
                cell.setDapiFociParams(fociNb, fociVol, fociInt, fociInt-cell.params.get("dapiBg")*fociVol/pixVol);
                break;
            case "GFP" :
                cell.setGfpFociParams(fociNb, fociVol, fociInt, fociInt-cell.params.get("gfpBg")*fociVol/pixVol);
                break;
        }
    }
    
    /**
     * Label object
     * @param popObj
     * @param img 
     * @param fontSize 
     */
    public void labelObject(Object3DInt obj, ImagePlus img, int fontSize) {
        BoundingBox bb = obj.getBoundingBox();
        int z = bb.zmin;
        int x = bb.xmin;
        int y = bb.ymin;
        img.setSlice(z+1);
        ImageProcessor ip = img.getProcessor();
        ip.setFont(new Font("SansSerif", Font.PLAIN, fontSize));
        ip.setColor(255);
        ip.drawString(Integer.toString((int)obj.getLabel()), x, y);
        img.updateAndDraw();
    }

    public void drawResults_wofoci(ImagePlus imgDAPI, ImagePlus imgPV, ArrayList<Cell> cells, String imageName, String outDir) {

        ImageHandler imgObj1 = ImageHandler.wrap(imgDAPI).createSameDimensions();
        ImageHandler imgObj2 = imgObj1.createSameDimensions();
        ImageHandler imgObj3 = imgObj1.createSameDimensions();
        ImageHandler imgObj4 = imgObj1.createSameDimensions(); // For membrane
    
        for (Cell cell : cells) {
            boolean isPV = false;
            if (cell.getPvCell() != null) {
                cell.getPvCell().drawObject(imgObj1);
                cell.getNucleus().drawObject(imgObj1, 0);
                isPV = true;
            }
            if (cell.getPnnCell() != null) {
                cell.getPnnCell().drawObject(imgObj2);
                cell.getNucleus().drawObject(imgObj2, 0);
            }
            cell.getNucleus().drawObject(imgObj3);
    
            // Draw membrane
            if (cell.getMembrane() != null) {
                cell.getMembrane().drawObject(imgObj4);
            }
    
            if (isPV)
                labelObject(cell.getPvCell(), imgObj3.getImagePlus(), 50);
            else
                labelObject(cell.getPnnCell(), imgObj3.getImagePlus(), 50);
        }
    
        ImagePlus[] imgColors1 = {imgObj1.getImagePlus(), imgObj2.getImagePlus(), imgObj3.getImagePlus(), imgObj4.getImagePlus(), imgPV};
        ImagePlus imgObjects1 = new RGBStackMerge().mergeHyperstacks(imgColors1, true);
        imgObjects1.setCalibration(imgDAPI.getCalibration());
        FileSaver ImgObjectsFile1 = new FileSaver(imgObjects1);
        ImgObjectsFile1.saveAsTiff(outDir + imageName + "_cells.tif");
        flush_close(imgObjects1);
    
        imgObj1.closeImagePlus();
        imgObj2.closeImagePlus();
        imgObj3.closeImagePlus();
        imgObj4.closeImagePlus();
    }

    /**
     * Save detected cells and foci in image
     */
    public void drawResults(ImagePlus imgDAPI, ImagePlus imgPV, ArrayList<Cell> cells, Objects3DIntPopulation gfpFociPop, 
        Objects3DIntPopulation dapiFociPop, String imageName, String outDir) {

    ImageHandler imgObj1 = ImageHandler.wrap(imgDAPI).createSameDimensions();
    ImageHandler imgObj2 = imgObj1.createSameDimensions();
    ImageHandler imgObj3 = imgObj1.createSameDimensions();
    ImageHandler imgObj4 = imgObj1.createSameDimensions(); // For membrane

    for (Cell cell : cells) {
        boolean isPV = false;
        if (cell.getPvCell() != null) {
            cell.getPvCell().drawObject(imgObj1);
            isPV = true;
        }
        if (cell.getPnnCell() != null) {
            cell.getPnnCell().drawObject(imgObj2);
        }
        cell.getNucleus().drawObject(imgObj3);

        // Draw membrane
        if (cell.getMembrane() != null) {
            cell.getMembrane().drawObject(imgObj4);
        }

        if (isPV)
            labelObject(cell.getPvCell(), imgObj3.getImagePlus(), 50);
        else
            labelObject(cell.getPnnCell(), imgObj3.getImagePlus(), 50);
    }

    ImagePlus[] imgColors1 = {imgObj1.getImagePlus(), imgObj2.getImagePlus(), imgObj3.getImagePlus(), imgObj4.getImagePlus(), imgPV};
    ImagePlus imgObjects1 = new RGBStackMerge().mergeHyperstacks(imgColors1, true);
    imgObjects1.setCalibration(imgDAPI.getCalibration());
    FileSaver ImgObjectsFile1 = new FileSaver(imgObjects1);
    ImgObjectsFile1.saveAsTiff(outDir + imageName + "_cells.tif");
    flush_close(imgObjects1);

    ImageHandler imgObj5 = imgObj1.createSameDimensions();
    ImageHandler imgObj6 = imgObj1.createSameDimensions();
    for (Object3DInt dot : gfpFociPop.getObjects3DInt()) {
        dot.drawObject(imgObj5, 255);
        dot.drawObject(imgObj3, 0);
    }
    for (Object3DInt dot : dapiFociPop.getObjects3DInt()) {
        dot.drawObject(imgObj6, 255);
        dot.drawObject(imgObj3, 0);
    }

    ImagePlus[] imgColors2 = {imgObj5.getImagePlus(), imgObj4.getImagePlus(), imgObj3.getImagePlus(), imgDAPI};
    ImagePlus imgObjects2 = new RGBStackMerge().mergeHyperstacks(imgColors2, false);
    imgObjects2.setCalibration(imgDAPI.getCalibration());
    FileSaver ImgObjectsFile2 = new FileSaver(imgObjects2);
    ImgObjectsFile2.saveAsTiff(outDir + imageName + "_foci.tif");
    flush_close(imgObjects2);

    imgObj1.closeImagePlus();
    imgObj2.closeImagePlus();
    imgObj3.closeImagePlus();
    imgObj4.closeImagePlus();
    imgObj5.closeImagePlus();
    imgObj6.closeImagePlus();
}

    public Objects3DIntPopulation stardistFociInMbPop(ImagePlus img, ArrayList<Cell> cells, String fociType, boolean resize) throws IOException{
        float fociIndex = 1;
        double resizeFactor = 0.5;
        Objects3DIntPopulation allFociPop = new Objects3DIntPopulation();
        for (Cell cell: cells) {
            Object3DInt mb = cell.getMembrane();
            if (mb == null) {
                continue; // Skip cells without a membrane
            }
            
            // Crop image around nucleus
            BoundingBox box = mb.getBoundingBox();
            Roi roiBox = new Roi(box.xmin, box.ymin, box.xmax-box.xmin, box.ymax-box.ymin);
            img.setRoi(roiBox);
            img.updateAndDraw();
            ImagePlus imgNuc = new Duplicator().run(img, box.zmin+1, box.zmax+1);
            imgNuc.deleteRoi();
            imgNuc.updateAndDraw();
            
            // Median filter, downscaling and gaussian filter
            ImagePlus imgM = median_filter(imgNuc, 1, 1);
            ImagePlus imgS = (resize) ? imgM.resize((int)(resizeFactor*imgNuc.getWidth()), (int)(resizeFactor*imgNuc.getHeight()), 1, "average") : imgM.duplicate();
            flush_close(imgM);
            ImagePlus imgG = gaussian_filter(imgS,1, 1);
            flush_close(imgS);

            // StarDist
            File starDistModelFile = new File(modelsPath+File.separator+stardistFociModel);
            StarDist2D star = new StarDist2D(syncObject, starDistModelFile);
            System.out.println(starDistModelFile +"path for stardist in membrane ");
            star.loadInput(imgG);
            star.setParams(stardistPercentileBottom, stardistPercentileTop, stardistFociProbThresh, stardistFociOverlayThresh, stardistOutput);
            star.run();
            flush_close(imgG);

            // Label foci in 3D
            ImagePlus imgLabels = star.associateLabels();
            if (resize) imgLabels = imgLabels.resize(imgNuc.getWidth(), imgNuc.getHeight(), 1, "none");
            
            flush_close(imgNuc);
            imgLabels.setCalibration(cal);
            Objects3DIntPopulation fociPop = new Objects3DIntPopulation(ImageHandler.wrap(imgLabels));
            popFilterSize(fociPop, minFociVol, maxFociVol);
            flush_close(imgLabels);
            
            // Find foci in nucleus
            fociPop.translateObjects(box.xmin, box.ymin, box.zmin);
            Objects3DIntPopulation fociColocPop = findFociInCell(mb, fociPop);
            System.out.println(fociColocPop.getNbObjects() + " " + fociType + " foci found in membrane " + mb.getLabel());
            
            for (Object3DInt foci: fociColocPop.getObjects3DInt()) {
                foci.setLabel(fociIndex);
                fociIndex++;
                foci.setType(cell.params.get("label").intValue());
                allFociPop.addObject(foci);
            }
            
            writeFociParameters(cell, fociColocPop, img, fociType);
        }
        return(allFociPop);
    }
} 
