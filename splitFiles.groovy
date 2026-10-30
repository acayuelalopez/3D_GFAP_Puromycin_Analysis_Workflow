import ij.IJ
import ij.ImagePlus
import ij.ImageStack
import ij.WindowManager
import ij.measure.Calibration
import ij.measure.ResultsTable
import ij.plugin.ChannelSplitter
import ij.plugin.Duplicator
import ij.plugin.RGBStackMerge
import ij.plugin.ZProjector
import inra.ijpb.label.LabelImages
import inra.ijpb.measure.region3d.RegionAnalyzer3D
import inra.ijpb.morphology.Strel
import loci.plugins.BF
import loci.plugins.in.ImporterOptions
import org.apache.commons.math3.stat.inference.KolmogorovSmirnovTest
import java.util.stream.Stream;

// INPUT UI
//
//#@File(label = "Input LIF File Directory", style = "directory") inputFiles
//#@File(label = "Output directory", style = "directory") outputDir
//#@File(label = "Green model", style = "file") greenModel
//#@File(label = "Red model", style = "file") redModel
//#@Integer(label = "Reference Channel", value = 1) refIndex
//#@Integer(label = "Target Channel", value = 2) targetIndex
def inputFiles = new File("/mnt/imgserver/IA/Projects/2025/2025_12_16_jbenzal/jbenzal/sub")
def outputDir = new File("/mnt/imgserver/IA/Projects/2025/2025_12_16_jbenzal/jbenzal/images_sub")

// IDE
//
//
//def headless = true;
//new ImageJ().setVisible(true);

IJ.log("-Parameters selected: ")
IJ.log("    -inputFileDir Ref: " + inputFiles)
IJ.log("    -outputDir: " + outputDir)
//IJ.log("    -Green Model: "+greenModel)
//IJ.log("    -Red Model: "+redModel)


IJ.log("                                                           ");
/** Get files (images) from input directory */
def listOfFiles = inputFiles.listFiles(); ;

for (def i = 0; i < listOfFiles.length; i++) {
/** Importer options for .lif file */
    def options = new ImporterOptions();
    options.setId(inputFiles.getAbsolutePath() + File.separator + listOfFiles[i].getName());
    options.setSplitChannels(false);
    options.setSplitTimepoints(false);
    options.setSplitFocalPlanes(false);
    options.setAutoscale(true);
    options.setStackFormat(ImporterOptions.VIEW_HYPERSTACK);
    options.setStackOrder(ImporterOptions.ORDER_XYCZT);
    options.setColorMode(ImporterOptions.COLOR_MODE_COMPOSITE);
    options.setCrop(false);
    options.setOpenAllSeries(true);
    def imps = BF.openImagePlus(options);

    for (def j = 0; j < imps.length; j++) {
        //if (!imps[j].getTitle().contains("JC ") && imps[j].getNSlices() > 1) {
            /** Get image serie per lif */
            def imp = imps[j]
            IJ.saveAs(imp, "Tiff", outputDir.getAbsolutePath() +File.separator+ imp.getTitle().replaceAll("/", ""))
            //def channels = ChannelSplitter.split(imp)
            //def impCP = RGBStackMerge.mergeChannels(new ImagePlus[]{channels[0], channels[1], channels[3]},false)
            //IJ.saveAs(impCP, "Tiff", outputDir.getAbsolutePath() + File.separator +"CP"+File.separator+ imp.getTitle().replaceAll("/", ""))
        }


    //}
}