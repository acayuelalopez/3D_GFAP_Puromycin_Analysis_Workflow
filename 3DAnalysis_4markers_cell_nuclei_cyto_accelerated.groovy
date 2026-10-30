import ij.IJ
import ij.ImagePlus
import ij.measure.ResultsTable
import ij.plugin.ChannelSplitter

import mcib3d.geom.Object3D
import mcib3d.geom.Object3DVoxels
import mcib3d.geom.Objects3DPopulation
import mcib3d.geom.Voxel3D
import mcib3d.image3d.ImageInt

import java.io.File
import java.util.ArrayList
import java.util.LinkedHashMap
import java.util.LinkedList
import java.util.List
import java.util.Locale
import java.util.Map


// =====================================================
// DIRECTORIOS
// =====================================================

final File inputFilesNuclei = new File(
        "/mnt/imgserver/IA/Projects/2025/2025_12_12_npriego/output/nuclei"
)

final File inputFilesRaw = new File(
        "/mnt/imgserver/IA/Projects/2025/2025_12_12_npriego/output/images"
)

final File outputDir = new File(
        "/mnt/imgserver/IA/Projects/2025/2025_12_12_npriego/output/csv"
)

if (!inputFilesNuclei.exists()) {
    IJ.log("ERROR: no existe el directorio de máscaras:")
    IJ.log(inputFilesNuclei.getAbsolutePath())
    return
}

if (!inputFilesRaw.exists()) {
    IJ.log("ERROR: no existe el directorio de imágenes:")
    IJ.log(inputFilesRaw.getAbsolutePath())
    return
}

if (!outputDir.exists() && !outputDir.mkdirs()) {
    IJ.log("ERROR: no se pudo crear el directorio de salida:")
    IJ.log(outputDir.getAbsolutePath())
    return
}


// =====================================================
// CONFIGURACIÓN DE CANALES
//
// Los índices empiezan en cero:
//
// Canal 1 de Fiji = índice 0
// Canal 2 de Fiji = índice 1
// Canal 3 de Fiji = índice 2
// Canal 4 de Fiji = índice 3
// =====================================================

final int DAPI_CHANNEL = 3
final int PUROMYCIN_CHANNEL = 1
final int GFAP_CHANNEL = 0


// =====================================================
// PARÁMETROS
// =====================================================

/*
 * Radio de dilatación expresado en píxeles XY.
 *
 * El radio Z se calcula utilizando la calibración física.
 */
final float CYTO_DILATION_XY_PX = 5.0f

final double GFAP_THRESHOLD_K = 0.5d
final double PURO_THRESHOLD_K = 0.5d


/*
 * false:
 *     Conserva la lógica original:
 *
 *     - Umbral calculado con medias celulares.
 *     - Clasificación realizada con medias nucleares.
 *
 * true:
 *     - Umbral calculado con medias celulares.
 *     - Clasificación realizada con medias celulares.
 */
final boolean CLASSIFY_USING_CELL_MEAN = false


/*
 * Si es true, ningún citoplasma virtual podrá contener
 * vóxeles pertenecientes a núcleos segmentados.
 */
final boolean EXCLUDE_ALL_NUCLEI_FROM_CYTO = true


/*
 * Si es true, guarda los CSV después de cada imagen
 * procesada correctamente.
 */
final boolean SAVE_AFTER_EACH_IMAGE = true


// =====================================================
// ARCHIVOS DE SALIDA
// =====================================================

final File perNucleusOutput = new File(
        outputDir,
        "3DAnalysis_perNucleus_GFAP_Puromycin_virtualCyto3D.csv"
)

final File perImageOutput = new File(
        outputDir,
        "3DAnalysis_perImage_GFAP_Puromycin_virtualCyto3D.csv"
)


// =====================================================
// TABLAS DE RESULTADOS
// =====================================================

final ResultsTable rtPerNucleus = new ResultsTable()
final ResultsTable rtPerImage = new ResultsTable()


// =====================================================
// OBTENER ARCHIVOS TIFF
// =====================================================

File[] allFiles = inputFilesRaw.listFiles()

if (allFiles == null) {
    IJ.log("ERROR: no se pudieron listar los archivos de:")
    IJ.log(inputFilesRaw.getAbsolutePath())
    return
}

List<File> validFiles = allFiles.findAll { File file ->

    if (file == null || !file.isFile()) {
        return false
    }

    String lowerName = file.getName().toLowerCase(Locale.ROOT)

    return !lowerName.contains(".ds_store") &&
            !lowerName.contains("_seg.npy") &&
            (
                    lowerName.endsWith(".tif") ||
                            lowerName.endsWith(".tiff")
            )

}.sort { File first, File second ->

    first.getName() <=> second.getName()
}

if (validFiles.isEmpty()) {
    IJ.log("ERROR: no se encontraron imágenes TIFF en:")
    IJ.log(inputFilesRaw.getAbsolutePath())
    return
}


// =====================================================
// INFORMACIÓN INICIAL
// =====================================================

long completeAnalysisStart = System.currentTimeMillis()

int processedImages = 0
int failedImages = 0

IJ.log("")
IJ.log("================================================")
IJ.log("ANÁLISIS 3D OPTIMIZADO")
IJ.log("================================================")
IJ.log("Número de imágenes: " + validFiles.size())
IJ.log("Radio XY: " + CYTO_DILATION_XY_PX + " píxeles")
IJ.log("Dilatación: 3D calibrada")
IJ.log(
        "Clasificación con media celular: " +
                CLASSIFY_USING_CELL_MEAN
)
IJ.log(
        "Excluir núcleos del citoplasma: " +
                EXCLUDE_ALL_NUCLEI_FROM_CYTO
)
IJ.log("================================================")


// =====================================================
// BUCLE PRINCIPAL
// =====================================================

for (int imageIndex = 0;
     imageIndex < validFiles.size();
     imageIndex++) {

    File rawFile = validFiles.get(imageIndex)

    long imageStart = System.currentTimeMillis()

    ImagePlus rawImage = null
    ImagePlus nucleiLabels = null
    ImagePlus[] channels = null

    boolean imageCompletedSuccessfully = false

    IJ.log("")
    IJ.log("-----------------------------------------------")
    IJ.log(
            "Procesando imagen " +
                    (imageIndex + 1) +
                    "/" +
                    validFiles.size()
    )
    IJ.log(rawFile.getName())

    try {

        // =================================================
        // LOCALIZAR MÁSCARA NUCLEAR
        // =================================================

        String nucleiMaskName =
                createMaskFileName(rawFile.getName())

        File nucleiMaskFile = new File(
                inputFilesNuclei,
                nucleiMaskName
        )

        if (!nucleiMaskFile.exists()) {
            IJ.log("AVISO: no existe la máscara:")
            IJ.log(nucleiMaskFile.getAbsolutePath())
            failedImages++
            continue
        }


        // =================================================
        // ABRIR IMÁGENES
        // =================================================

        nucleiLabels = IJ.openImage(
                nucleiMaskFile.getAbsolutePath()
        )

        if (nucleiLabels == null ||
                nucleiLabels.getProcessor() == null) {

            IJ.log(
                    "AVISO: no se pudo abrir la máscara nuclear."
            )

            failedImages++
            continue
        }

        rawImage = IJ.openImage(
                rawFile.getAbsolutePath()
        )

        if (rawImage == null ||
                rawImage.getProcessor() == null) {

            IJ.log(
                    "AVISO: no se pudo abrir la imagen raw."
            )

            failedImages++
            continue
        }


        // =================================================
        // SEPARAR CANALES
        // =================================================

        channels = ChannelSplitter.split(rawImage)

        if (channels == null || channels.length == 0) {
            IJ.log(
                    "AVISO: no se pudieron separar los canales."
            )

            failedImages++
            continue
        }

        int maximumRequestedChannel = DAPI_CHANNEL

        if (PUROMYCIN_CHANNEL > maximumRequestedChannel) {
            maximumRequestedChannel = PUROMYCIN_CHANNEL
        }

        if (GFAP_CHANNEL > maximumRequestedChannel) {
            maximumRequestedChannel = GFAP_CHANNEL
        }

        if (maximumRequestedChannel >= channels.length) {

            IJ.log(
                    "ERROR: índice de canal fuera de rango."
            )

            IJ.log(
                    "Número de canales disponibles: " +
                            channels.length
            )

            IJ.log(
                    "Índice máximo solicitado: " +
                            maximumRequestedChannel
            )

            failedImages++
            continue
        }

        ImagePlus channelDapi =
                channels[DAPI_CHANNEL]

        ImagePlus channelPuromycin =
                channels[PUROMYCIN_CHANNEL]

        ImagePlus channelGfap =
                channels[GFAP_CHANNEL]


        // =================================================
        // COMPROBAR CANALES
        // =================================================

        if (channelDapi == null ||
                channelPuromycin == null ||
                channelGfap == null) {

            IJ.log(
                    "ERROR: uno o varios canales son nulos."
            )

            failedImages++
            continue
        }


        // =================================================
        // COMPROBAR DIMENSIONES
        // =================================================

        if (!sameDimensions(
                nucleiLabels,
                channelGfap
        )) {

            IJ.log(
                    "ERROR: dimensiones diferentes entre máscara y GFAP."
            )

            IJ.log(
                    "Máscara: " +
                            dimensionsAsString(nucleiLabels)
            )

            IJ.log(
                    "GFAP: " +
                            dimensionsAsString(channelGfap)
            )

            failedImages++
            continue
        }

        if (!sameDimensions(
                nucleiLabels,
                channelPuromycin
        )) {

            IJ.log(
                    "ERROR: dimensiones diferentes entre máscara y puromicina."
            )

            IJ.log(
                    "Máscara: " +
                            dimensionsAsString(nucleiLabels)
            )

            IJ.log(
                    "Puromicina: " +
                            dimensionsAsString(channelPuromycin)
            )

            failedImages++
            continue
        }


        // =================================================
        // CONFIGURAR CALIBRACIÓN
        // =================================================

        nucleiLabels.setCalibration(
                channelDapi.getCalibration().copy()
        )

        double pixelWidth =
                nucleiLabels.getCalibration().pixelWidth as double

        double pixelHeight =
                nucleiLabels.getCalibration().pixelHeight as double

        double pixelDepth =
                nucleiLabels.getCalibration().pixelDepth as double

        if (pixelWidth <= 0.0d) {
            pixelWidth = 1.0d
        }

        if (pixelHeight <= 0.0d) {
            pixelHeight = pixelWidth
        }

        if (pixelDepth <= 0.0d) {
            pixelDepth = pixelWidth
        }

        if (Math.abs(pixelWidth - pixelHeight) > 1.0E-9d) {

            IJ.log(
                    "AVISO: pixelWidth y pixelHeight son distintos."
            )

            IJ.log(
                    "pixelWidth = " +
                            pixelWidth +
                            ", pixelHeight = " +
                            pixelHeight
            )
        }

        /*
         * Distancia física correspondiente al radio XY.
         */
        double dilationPhysicalDistance =
                ((double) CYTO_DILATION_XY_PX) *
                        pixelWidth

        /*
         * Radio de dilatación en coordenadas de vóxel.
         *
         * No se fuerza un mínimo de un slice Z, porque eso
         * alteraría la distancia física solicitada.
         */
        float dilationRadiusX =
                CYTO_DILATION_XY_PX

        float dilationRadiusY =
                (float) (
                        dilationPhysicalDistance /
                                pixelHeight
                )

        float dilationRadiusZ =
                (float) (
                        dilationPhysicalDistance /
                                pixelDepth
                )

        IJ.log(
                "Calibración XYZ: " +
                        pixelWidth +
                        " x " +
                        pixelHeight +
                        " x " +
                        pixelDepth
        )

        IJ.log(
                "Radio 3D XYZ en vóxeles: " +
                        dilationRadiusX +
                        " x " +
                        dilationRadiusY +
                        " x " +
                        dilationRadiusZ
        )

        IJ.log(
                "Distancia física de dilatación: " +
                        dilationPhysicalDistance +
                        " " +
                        nucleiLabels
                                .getCalibration()
                                .getUnit()
        )

        if (dilationRadiusZ < 1.0f) {

            IJ.log(
                    "NOTA: el radio físico es menor que un slice Z."
            )

            IJ.log(
                    "La dilatación está configurada en 3D, pero puede " +
                            "no alcanzar el slice Z vecino."
            )
        }


        // =================================================
        // CREAR OBJETOS MCIB3D
        // =================================================

        ImageInt imageNuclei =
                ImageInt.wrap(nucleiLabels)

        ImageInt signalGfap =
                ImageInt.wrap(channelGfap)

        ImageInt signalPuromycin =
                ImageInt.wrap(channelPuromycin)

        Objects3DPopulation nucleiPopulation =
                new Objects3DPopulation(imageNuclei)

        int numberOfNuclei =
                nucleiPopulation.getNbObjects()

        if (numberOfNuclei == 0) {
            IJ.log(
                    "AVISO: la máscara no contiene núcleos."
            )

            failedImages++
            continue
        }

        IJ.log(
                "Número de núcleos: " +
                        numberOfNuclei
        )


        // =================================================
        // UNA ÚNICA PASADA COSTOSA DE MEDICIÓN
        // =================================================

        List<Map<String, Object>> measurements =
                new ArrayList<Map<String, Object>>(
                        numberOfNuclei
                )

        double gfapCellMeanSum = 0.0d
        double gfapCellMeanSumSquares = 0.0d

        double puroCellMeanSum = 0.0d
        double puroCellMeanSumSquares = 0.0d


        /*
         * CORRECCIÓN DEL ERROR:
         *
         * No usar:
         *
         * Math.max(1, numberOfNuclei / 10)
         *
         * porque Groovy puede convertir la división en BigDecimal.
         */
        int progressStep =
                numberOfNuclei.intdiv(10)

        if (progressStep < 1) {
            progressStep = 1
        }


        for (int nucleusIndex = 0;
             nucleusIndex < numberOfNuclei;
             nucleusIndex++) {

            if (nucleusIndex % progressStep == 0) {

                IJ.showProgress(
                        nucleusIndex,
                        numberOfNuclei
                )

                IJ.showStatus(
                        "Midiendo núcleos de " +
                                rawFile.getName() +
                                ": " +
                                nucleusIndex +
                                "/" +
                                numberOfNuclei
                )
            }

            Object3D nucleus =
                    nucleiPopulation.getObject(nucleusIndex)


            // ---------------------------------------------
            // MEDIDAS DEL NÚCLEO
            // ---------------------------------------------

            double nucleusVolume =
                    nucleus.getVolumePixels() as double

            double gfapSumNucleus =
                    nucleus.getIntegratedDensity(
                            signalGfap
                    ) as double

            double puroSumNucleus =
                    nucleus.getIntegratedDensity(
                            signalPuromycin
                    ) as double

            double gfapMeanNucleus =
                    nucleus.getPixMeanValue(
                            signalGfap
                    ) as double

            double puroMeanNucleus =
                    nucleus.getPixMeanValue(
                            signalPuromycin
                    ) as double


            // ---------------------------------------------
            // DILATACIÓN 3D NATIVA DE MCIB3D
            // ---------------------------------------------

            Object3DVoxels dilatedObject =
                    nucleus.getDilatedObject(
                            dilationRadiusX,
                            dilationRadiusY,
                            dilationRadiusZ
                    )

            Object3DVoxels virtualCytoObject = null

            if (dilatedObject != null) {

                /*
                 * Citoplasma virtual:
                 *
                 * objeto dilatado - núcleo original
                 */
                dilatedObject.substractObject(nucleus)

                /*
                 * Elimina los vóxeles fuera de la imagen y,
                 * opcionalmente, todos los vóxeles nucleares.
                 */
                virtualCytoObject =
                        cleanVirtualCytoObject(
                                dilatedObject,
                                imageNuclei,
                                nucleiLabels.getWidth(),
                                nucleiLabels.getHeight(),
                                nucleiLabels.getNSlices(),
                                EXCLUDE_ALL_NUCLEI_FROM_CYTO,
                                pixelWidth,
                                pixelDepth,
                                nucleiLabels
                                        .getCalibration()
                                        .getUnit()
                        )
            }


            // ---------------------------------------------
            // MEDIDAS DEL CITOPLASMA VIRTUAL
            // ---------------------------------------------

            double virtualCytoVolume = 0.0d

            double gfapSumVirtualCyto = 0.0d
            double gfapMeanVirtualCyto = 0.0d

            double puroSumVirtualCyto = 0.0d
            double puroMeanVirtualCyto = 0.0d

            if (virtualCytoObject != null) {

                virtualCytoVolume =
                        virtualCytoObject
                                .getVolumePixels() as double

                if (virtualCytoVolume > 0.0d) {

                    gfapSumVirtualCyto =
                            virtualCytoObject
                                    .getIntegratedDensity(
                                            signalGfap
                                    ) as double

                    puroSumVirtualCyto =
                            virtualCytoObject
                                    .getIntegratedDensity(
                                            signalPuromycin
                                    ) as double

                    gfapMeanVirtualCyto =
                            gfapSumVirtualCyto /
                                    virtualCytoVolume

                    puroMeanVirtualCyto =
                            puroSumVirtualCyto /
                                    virtualCytoVolume
                }
            }


            // ---------------------------------------------
            // MEDIDAS DE LA CÉLULA
            // ---------------------------------------------

            /*
             * No se construye una máscara celular.
             *
             * El núcleo y el citoplasma son disjuntos:
             *
             * célula = núcleo + citoplasma
             */
            double cellVolume =
                    nucleusVolume +
                            virtualCytoVolume

            double gfapSumCell =
                    gfapSumNucleus +
                            gfapSumVirtualCyto

            double puroSumCell =
                    puroSumNucleus +
                            puroSumVirtualCyto

            double gfapMeanCell = 0.0d
            double puroMeanCell = 0.0d

            if (cellVolume > 0.0d) {

                gfapMeanCell =
                        gfapSumCell /
                                cellVolume

                puroMeanCell =
                        puroSumCell /
                                cellVolume
            }


            // ---------------------------------------------
            // ACUMULAR PARA LOS UMBRALES
            // ---------------------------------------------

            gfapCellMeanSum +=
                    gfapMeanCell

            gfapCellMeanSumSquares +=
                    gfapMeanCell *
                            gfapMeanCell

            puroCellMeanSum +=
                    puroMeanCell

            puroCellMeanSumSquares +=
                    puroMeanCell *
                            puroMeanCell


            // ---------------------------------------------
            // GUARDAR MEDIDAS
            // ---------------------------------------------

            Map<String, Object> measurement =
                    new LinkedHashMap<String, Object>()

            measurement.put(
                    "nucleus",
                    nucleus
            )

            measurement.put(
                    "nucleusVolume",
                    nucleusVolume
            )

            measurement.put(
                    "gfapMeanNucleus",
                    gfapMeanNucleus
            )

            measurement.put(
                    "gfapSumNucleus",
                    gfapSumNucleus
            )

            measurement.put(
                    "puroMeanNucleus",
                    puroMeanNucleus
            )

            measurement.put(
                    "puroSumNucleus",
                    puroSumNucleus
            )

            measurement.put(
                    "virtualCytoVolume",
                    virtualCytoVolume
            )

            measurement.put(
                    "gfapMeanVirtualCyto",
                    gfapMeanVirtualCyto
            )

            measurement.put(
                    "gfapSumVirtualCyto",
                    gfapSumVirtualCyto
            )

            measurement.put(
                    "puroMeanVirtualCyto",
                    puroMeanVirtualCyto
            )

            measurement.put(
                    "puroSumVirtualCyto",
                    puroSumVirtualCyto
            )

            measurement.put(
                    "cellVolume",
                    cellVolume
            )

            measurement.put(
                    "gfapMeanCell",
                    gfapMeanCell
            )

            measurement.put(
                    "gfapSumCell",
                    gfapSumCell
            )

            measurement.put(
                    "puroMeanCell",
                    puroMeanCell
            )

            measurement.put(
                    "puroSumCell",
                    puroSumCell
            )

            measurements.add(measurement)
        }

        IJ.showProgress(1.0d)


        // =================================================
        // CALCULAR UMBRALES
        // =================================================

        int measuredNuclei =
                measurements.size()

        if (measuredNuclei == 0) {
            IJ.log(
                    "AVISO: no se pudo medir ningún núcleo."
            )

            failedImages++
            continue
        }

        double measuredNucleiDouble =
                (double) measuredNuclei

        double gfapMeanAll =
                gfapCellMeanSum /
                        measuredNucleiDouble

        double puroMeanAll =
                puroCellMeanSum /
                        measuredNucleiDouble

        double gfapVariance =
                gfapCellMeanSumSquares /
                        measuredNucleiDouble -
                        gfapMeanAll *
                        gfapMeanAll

        double puroVariance =
                puroCellMeanSumSquares /
                        measuredNucleiDouble -
                        puroMeanAll *
                        puroMeanAll

        /*
         * No usamos Math.max para evitar ambigüedades
         * de tipos en versiones antiguas de Groovy.
         */
        if (gfapVariance < 0.0d) {
            gfapVariance = 0.0d
        }

        if (puroVariance < 0.0d) {
            puroVariance = 0.0d
        }

        double gfapStandardDeviation =
                Math.sqrt(gfapVariance)

        double puroStandardDeviation =
                Math.sqrt(puroVariance)

        double gfapThreshold =
                gfapMeanAll +
                        GFAP_THRESHOLD_K *
                        gfapStandardDeviation

        double puroThreshold =
                puroMeanAll +
                        PURO_THRESHOLD_K *
                        puroStandardDeviation

        IJ.log(
                "Umbral GFAP: " +
                        gfapThreshold
        )

        IJ.log(
                "Umbral puromicina: " +
                        puroThreshold
        )


        // =================================================
        // ACUMULADORES POR IMAGEN
        // =================================================

        int nNuclei = measuredNuclei
        int nGFAPpos = 0
        int nGFAPposPuropos = 0
        int nDapiGfapPos = 0
        int nDapiGfapPuroPos = 0

        List<Double> puroMeanGFAPposMasks =
                new ArrayList<Double>()

        List<Double> puroSumGFAPposMasks =
                new ArrayList<Double>()

        List<Double> gfapMeanDapiGfap =
                new ArrayList<Double>()

        List<Double> gfapSumDapiGfap =
                new ArrayList<Double>()

        List<Double> puroMeanDapiGfap =
                new ArrayList<Double>()

        List<Double> puroSumDapiGfap =
                new ArrayList<Double>()

        List<Double> gfapMeanDapiGfapPuro =
                new ArrayList<Double>()

        List<Double> gfapSumDapiGfapPuro =
                new ArrayList<Double>()

        List<Double> puroMeanDapiGfapPuro =
                new ArrayList<Double>()

        List<Double> puroSumDapiGfapPuro =
                new ArrayList<Double>()

        List<Double> gfapMeanVirtualCytoDapiGfap =
                new ArrayList<Double>()

        List<Double> gfapSumVirtualCytoDapiGfap =
                new ArrayList<Double>()

        List<Double> puroMeanVirtualCytoDapiGfap =
                new ArrayList<Double>()

        List<Double> puroSumVirtualCytoDapiGfap =
                new ArrayList<Double>()

        List<Double> gfapMeanCellDapiGfap =
                new ArrayList<Double>()

        List<Double> gfapSumCellDapiGfap =
                new ArrayList<Double>()

        List<Double> puroMeanCellDapiGfap =
                new ArrayList<Double>()

        List<Double> puroSumCellDapiGfap =
                new ArrayList<Double>()


        // =================================================
        // CLASIFICAR SIN VOLVER A CREAR MÁSCARAS
        // =================================================

        for (Map<String, Object> measurement :
                measurements) {

            Object3D nucleus =
                    measurement.get("nucleus") as Object3D

            double nucleusVolume =
                    numberValue(
                            measurement.get("nucleusVolume")
                    )

            double gfapMeanNucleus =
                    numberValue(
                            measurement.get("gfapMeanNucleus")
                    )

            double gfapSumNucleus =
                    numberValue(
                            measurement.get("gfapSumNucleus")
                    )

            double puroMeanNucleus =
                    numberValue(
                            measurement.get("puroMeanNucleus")
                    )

            double puroSumNucleus =
                    numberValue(
                            measurement.get("puroSumNucleus")
                    )

            double virtualCytoVolume =
                    numberValue(
                            measurement.get(
                                    "virtualCytoVolume"
                            )
                    )

            double gfapMeanVirtualCyto =
                    numberValue(
                            measurement.get(
                                    "gfapMeanVirtualCyto"
                            )
                    )

            double gfapSumVirtualCyto =
                    numberValue(
                            measurement.get(
                                    "gfapSumVirtualCyto"
                            )
                    )

            double puroMeanVirtualCyto =
                    numberValue(
                            measurement.get(
                                    "puroMeanVirtualCyto"
                            )
                    )

            double puroSumVirtualCyto =
                    numberValue(
                            measurement.get(
                                    "puroSumVirtualCyto"
                            )
                    )

            double cellVolume =
                    numberValue(
                            measurement.get("cellVolume")
                    )

            double gfapMeanCell =
                    numberValue(
                            measurement.get("gfapMeanCell")
                    )

            double gfapSumCell =
                    numberValue(
                            measurement.get("gfapSumCell")
                    )

            double puroMeanCell =
                    numberValue(
                            measurement.get("puroMeanCell")
                    )

            double puroSumCell =
                    numberValue(
                            measurement.get("puroSumCell")
                    )


            // ---------------------------------------------
            // CLASIFICACIÓN
            // ---------------------------------------------

            boolean isDapiPos = true

            double gfapClassificationValue

            double puroClassificationValue

            if (CLASSIFY_USING_CELL_MEAN) {

                gfapClassificationValue =
                        gfapMeanCell

                puroClassificationValue =
                        puroMeanCell

            } else {

                gfapClassificationValue =
                        gfapMeanNucleus

                puroClassificationValue =
                        puroMeanNucleus
            }

            boolean isGfapPos =
                    gfapClassificationValue >
                            gfapThreshold

            boolean isPuroPos =
                    puroClassificationValue >
                            puroThreshold

            boolean isDapiGfapPos =
                    isDapiPos &&
                            isGfapPos

            boolean isDapiGfapPuroPos =
                    isDapiPos &&
                            isGfapPos &&
                            isPuroPos


            // ---------------------------------------------
            // CONTADORES GFAP
            // ---------------------------------------------

            if (isGfapPos) {

                nGFAPpos++

                puroMeanGFAPposMasks.add(
                        puroMeanNucleus
                )

                puroSumGFAPposMasks.add(
                        puroSumNucleus
                )
            }

            if (isGfapPos && isPuroPos) {
                nGFAPposPuropos++
            }


            // ---------------------------------------------
            // DAPI+ GFAP+
            // ---------------------------------------------

            if (isDapiGfapPos) {

                nDapiGfapPos++

                gfapMeanDapiGfap.add(
                        gfapMeanNucleus
                )

                gfapSumDapiGfap.add(
                        gfapSumNucleus
                )

                puroMeanDapiGfap.add(
                        puroMeanNucleus
                )

                puroSumDapiGfap.add(
                        puroSumNucleus
                )

                gfapMeanCellDapiGfap.add(
                        gfapMeanCell
                )

                gfapSumCellDapiGfap.add(
                        gfapSumCell
                )

                puroMeanCellDapiGfap.add(
                        puroMeanCell
                )

                puroSumCellDapiGfap.add(
                        puroSumCell
                )

                gfapMeanVirtualCytoDapiGfap.add(
                        gfapMeanVirtualCyto
                )

                gfapSumVirtualCytoDapiGfap.add(
                        gfapSumVirtualCyto
                )

                puroMeanVirtualCytoDapiGfap.add(
                        puroMeanVirtualCyto
                )

                puroSumVirtualCytoDapiGfap.add(
                        puroSumVirtualCyto
                )
            }


            // ---------------------------------------------
            // DAPI+ GFAP+ PUROMYCIN+
            // ---------------------------------------------

            if (isDapiGfapPuroPos) {

                nDapiGfapPuroPos++

                gfapMeanDapiGfapPuro.add(
                        gfapMeanNucleus
                )

                gfapSumDapiGfapPuro.add(
                        gfapSumNucleus
                )

                puroMeanDapiGfapPuro.add(
                        puroMeanNucleus
                )

                puroSumDapiGfapPuro.add(
                        puroSumNucleus
                )
            }


            // ---------------------------------------------
            // TABLA POR NÚCLEO
            // ---------------------------------------------

            rtPerNucleus.incrementCounter()

            int row =
                    rtPerNucleus.getCounter() - 1

            rtPerNucleus.setValue(
                    "Image",
                    row,
                    rawFile.getName()
            )

            rtPerNucleus.setValue(
                    "Nucleus_ID",
                    row,
                    nucleus.getValue()
            )

            rtPerNucleus.setValue(
                    "Nucleus_CenterX",
                    row,
                    nucleus.getCenterX()
            )

            rtPerNucleus.setValue(
                    "Nucleus_CenterY",
                    row,
                    nucleus.getCenterY()
            )

            rtPerNucleus.setValue(
                    "Nucleus_CenterZ",
                    row,
                    nucleus.getCenterZ()
            )

            rtPerNucleus.setValue(
                    "Nucleus_ZRange",
                    row,
                    nucleus.getZmax() -
                            nucleus.getZmin()
            )

            rtPerNucleus.setValue(
                    "Nucleus_Volume_px",
                    row,
                    nucleusVolume
            )

            rtPerNucleus.setValue(
                    "GFAP_Mean_Nucleus",
                    row,
                    gfapMeanNucleus
            )

            rtPerNucleus.setValue(
                    "GFAP_Sum_Nucleus",
                    row,
                    gfapSumNucleus
            )

            rtPerNucleus.setValue(
                    "Puromycin_Mean_Nucleus",
                    row,
                    puroMeanNucleus
            )

            rtPerNucleus.setValue(
                    "Puromycin_Sum_Nucleus",
                    row,
                    puroSumNucleus
            )

            rtPerNucleus.setValue(
                    "Cell_Volume_px",
                    row,
                    cellVolume
            )

            rtPerNucleus.setValue(
                    "GFAP_Mean_Cell",
                    row,
                    gfapMeanCell
            )

            rtPerNucleus.setValue(
                    "GFAP_Sum_Cell",
                    row,
                    gfapSumCell
            )

            rtPerNucleus.setValue(
                    "Puromycin_Mean_Cell",
                    row,
                    puroMeanCell
            )

            rtPerNucleus.setValue(
                    "Puromycin_Sum_Cell",
                    row,
                    puroSumCell
            )

            rtPerNucleus.setValue(
                    "GFAP_Threshold",
                    row,
                    gfapThreshold
            )

            rtPerNucleus.setValue(
                    "Puromycin_Threshold",
                    row,
                    puroThreshold
            )

            rtPerNucleus.setValue(
                    "DAPI_Status",
                    row,
                    isDapiPos ? 1 : 0
            )

            rtPerNucleus.setValue(
                    "GFAP_Status",
                    row,
                    isGfapPos ? 1 : 0
            )

            rtPerNucleus.setValue(
                    "Puromycin_Status",
                    row,
                    isPuroPos ? 1 : 0
            )

            rtPerNucleus.setValue(
                    "GFAP+_Puromycin+",
                    row,
                    isGfapPos && isPuroPos ? 1 : 0
            )

            rtPerNucleus.setValue(
                    "DAPI+_GFAP+",
                    row,
                    isDapiGfapPos ? 1 : 0
            )

            rtPerNucleus.setValue(
                    "DAPI+_GFAP+_Puromycin+",
                    row,
                    isDapiGfapPuroPos ? 1 : 0
            )

            rtPerNucleus.setValue(
                    "VirtualCyto_Dilation_XY_px",
                    row,
                    CYTO_DILATION_XY_PX
            )

            rtPerNucleus.setValue(
                    "VirtualCyto_Dilation_Z_slices",
                    row,
                    dilationRadiusZ
            )

            rtPerNucleus.setValue(
                    "VirtualCyto_Dilation_Physical",
                    row,
                    dilationPhysicalDistance
            )

            rtPerNucleus.setValue(
                    "VirtualCyto_Volume_px",
                    row,
                    virtualCytoVolume
            )

            rtPerNucleus.setValue(
                    "GFAP_Mean_VirtualCyto",
                    row,
                    gfapMeanVirtualCyto
            )

            rtPerNucleus.setValue(
                    "GFAP_Sum_VirtualCyto",
                    row,
                    gfapSumVirtualCyto
            )

            rtPerNucleus.setValue(
                    "Puromycin_Mean_VirtualCyto",
                    row,
                    puroMeanVirtualCyto
            )

            rtPerNucleus.setValue(
                    "Puromycin_Sum_VirtualCyto",
                    row,
                    puroSumVirtualCyto
            )
        }


        // =================================================
        // TABLA POR IMAGEN
        // =================================================

        rtPerImage.incrementCounter()

        int rowImage =
                rtPerImage.getCounter() - 1

        rtPerImage.setValue(
                "Image",
                rowImage,
                rawFile.getName()
        )

        rtPerImage.setValue(
                "N_Nuclei",
                rowImage,
                nNuclei
        )

        rtPerImage.setValue(
                "GFAP_Threshold",
                rowImage,
                gfapThreshold
        )

        rtPerImage.setValue(
                "Puromycin_Threshold",
                rowImage,
                puroThreshold
        )

        rtPerImage.setValue(
                "VirtualCyto_Dilation_XY_px",
                rowImage,
                CYTO_DILATION_XY_PX
        )

        rtPerImage.setValue(
                "VirtualCyto_Dilation_Z_slices",
                rowImage,
                dilationRadiusZ
        )

        rtPerImage.setValue(
                "VirtualCyto_Dilation_Physical",
                rowImage,
                dilationPhysicalDistance
        )


        // ---------------------------------------------
        // GFAP+
        // ---------------------------------------------

        rtPerImage.setValue(
                "N_GFAP+",
                rowImage,
                nGFAPpos
        )

        rtPerImage.setValue(
                "N_GFAP+_Puromycin+",
                rowImage,
                nGFAPposPuropos
        )

        rtPerImage.setValue(
                "Puromycin_Mean_GFAP+Masks_Avg",
                rowImage,
                meanValues(puroMeanGFAPposMasks)
        )

        rtPerImage.setValue(
                "Puromycin_Sum_GFAP+Masks_Sum",
                rowImage,
                sumValues(puroSumGFAPposMasks)
        )


        // ---------------------------------------------
        // DAPI+ GFAP+
        // ---------------------------------------------

        rtPerImage.setValue(
                "N_DAPI+_GFAP+",
                rowImage,
                nDapiGfapPos
        )

        rtPerImage.setValue(
                "GFAP_Mean_DAPI+_GFAP+_Nuclei_Avg",
                rowImage,
                meanValues(gfapMeanDapiGfap)
        )

        rtPerImage.setValue(
                "GFAP_Sum_DAPI+_GFAP+_Nuclei_Sum",
                rowImage,
                sumValues(gfapSumDapiGfap)
        )

        rtPerImage.setValue(
                "Puromycin_Mean_DAPI+_GFAP+_Nuclei_Avg",
                rowImage,
                meanValues(puroMeanDapiGfap)
        )

        rtPerImage.setValue(
                "Puromycin_Sum_DAPI+_GFAP+_Nuclei_Sum",
                rowImage,
                sumValues(puroSumDapiGfap)
        )


        // ---------------------------------------------
        // DAPI+ GFAP+ PUROMYCIN+
        // ---------------------------------------------

        rtPerImage.setValue(
                "N_DAPI+_GFAP+_Puromycin+",
                rowImage,
                nDapiGfapPuroPos
        )

        rtPerImage.setValue(
                "GFAP_Mean_DAPI+_GFAP+_Puromycin+_Nuclei_Avg",
                rowImage,
                meanValues(gfapMeanDapiGfapPuro)
        )

        rtPerImage.setValue(
                "GFAP_Sum_DAPI+_GFAP+_Puromycin+_Nuclei_Sum",
                rowImage,
                sumValues(gfapSumDapiGfapPuro)
        )

        rtPerImage.setValue(
                "Puromycin_Mean_DAPI+_GFAP+_Puromycin+_Nuclei_Avg",
                rowImage,
                meanValues(puroMeanDapiGfapPuro)
        )

        rtPerImage.setValue(
                "Puromycin_Sum_DAPI+_GFAP+_Puromycin+_Nuclei_Sum",
                rowImage,
                sumValues(puroSumDapiGfapPuro)
        )


        // ---------------------------------------------
        // CÉLULA COMPLETA DAPI+ GFAP+
        // ---------------------------------------------

        rtPerImage.setValue(
                "GFAP_Mean_Cell_DAPI+_GFAP+_Avg",
                rowImage,
                meanValues(gfapMeanCellDapiGfap)
        )

        rtPerImage.setValue(
                "GFAP_Sum_Cell_DAPI+_GFAP+_Sum",
                rowImage,
                sumValues(gfapSumCellDapiGfap)
        )

        rtPerImage.setValue(
                "Puromycin_Mean_Cell_DAPI+_GFAP+_Avg",
                rowImage,
                meanValues(puroMeanCellDapiGfap)
        )

        rtPerImage.setValue(
                "Puromycin_Sum_Cell_DAPI+_GFAP+_Sum",
                rowImage,
                sumValues(puroSumCellDapiGfap)
        )


        // ---------------------------------------------
        // CITOPLASMA VIRTUAL DAPI+ GFAP+
        // ---------------------------------------------

        rtPerImage.setValue(
                "GFAP_Mean_VirtualCyto_DAPI+_GFAP+_Avg",
                rowImage,
                meanValues(
                        gfapMeanVirtualCytoDapiGfap
                )
        )

        rtPerImage.setValue(
                "GFAP_Sum_VirtualCyto_DAPI+_GFAP+_Sum",
                rowImage,
                sumValues(
                        gfapSumVirtualCytoDapiGfap
                )
        )

        rtPerImage.setValue(
                "Puromycin_Mean_VirtualCyto_DAPI+_GFAP+_Avg",
                rowImage,
                meanValues(
                        puroMeanVirtualCytoDapiGfap
                )
        )

        rtPerImage.setValue(
                "Puromycin_Sum_VirtualCyto_DAPI+_GFAP+_Sum",
                rowImage,
                sumValues(
                        puroSumVirtualCytoDapiGfap
                )
        )


        // =================================================
        // IMAGEN COMPLETADA
        // =================================================

        imageCompletedSuccessfully = true
        processedImages++


        // =================================================
        // GUARDADO PROGRESIVO
        // =================================================

        if (SAVE_AFTER_EACH_IMAGE) {

            rtPerNucleus.saveAs(
                    perNucleusOutput.getAbsolutePath()
            )

            rtPerImage.saveAs(
                    perImageOutput.getAbsolutePath()
            )
        }


        // =================================================
        // TIEMPO POR IMAGEN
        // =================================================

        double imageSeconds =
                (
                        System.currentTimeMillis() -
                                imageStart
                ) / 1000.0d

        IJ.log(
                "Imagen terminada correctamente en " +
                        String.format(
                                Locale.US,
                                "%.2f",
                                imageSeconds
                        ) +
                        " segundos."
        )

        IJ.log(
                "Núcleos: " +
                        nNuclei +
                        " | GFAP+: " +
                        nGFAPpos +
                        " | GFAP+ Puromycin+: " +
                        nGFAPposPuropos
        )

    } catch (Throwable exception) {

        if (!imageCompletedSuccessfully) {
            failedImages++
        }

        IJ.log("")
        IJ.log(
                "ERROR procesando: " +
                        rawFile.getName()
        )

        IJ.log(exception.toString())

        for (StackTraceElement element :
                exception.getStackTrace()) {

            IJ.log(
                    "    " +
                            element.toString()
            )
        }

    } finally {

        // =================================================
        // LIBERAR MEMORIA ENTRE IMÁGENES
        // =================================================

        if (channels != null) {

            for (ImagePlus channel : channels) {

                if (channel != null) {
                    channel.close()
                }
            }
        }

        if (rawImage != null) {
            rawImage.close()
        }

        if (nucleiLabels != null) {
            nucleiLabels.close()
        }

        IJ.showProgress(1.0d)
    }
}


// =====================================================
// NO GUARDAR TABLAS VACÍAS
// =====================================================

if (rtPerNucleus.getCounter() == 0) {

    IJ.log("")
    IJ.log("================================================")
    IJ.log("ERROR: NO SE GENERARON RESULTADOS")
    IJ.log("================================================")
    IJ.log("No se guardarán CSV vacíos.")
    IJ.log("Imágenes procesadas correctamente: " + processedImages)
    IJ.log("Imágenes con error: " + failedImages)
    IJ.log("Revisa los errores anteriores.")
    IJ.log("================================================")

    return
}


// =====================================================
// GUARDADO FINAL
// =====================================================

try {

    rtPerNucleus.saveAs(
            perNucleusOutput.getAbsolutePath()
    )

    rtPerImage.saveAs(
            perImageOutput.getAbsolutePath()
    )

    double completeSeconds =
            (
                    System.currentTimeMillis() -
                            completeAnalysisStart
            ) / 1000.0d

    IJ.log("")
    IJ.log("================================================")
    IJ.log("ANÁLISIS TERMINADO")
    IJ.log("================================================")

    IJ.log(
            "Imágenes procesadas correctamente: " +
                    processedImages
    )

    IJ.log(
            "Imágenes con error: " +
                    failedImages
    )

    IJ.log(
            "Filas por núcleo: " +
                    rtPerNucleus.getCounter()
    )

    IJ.log(
            "Filas por imagen: " +
                    rtPerImage.getCounter()
    )

    IJ.log(
            "Tiempo total: " +
                    String.format(
                            Locale.US,
                            "%.2f",
                            completeSeconds
                    ) +
                    " segundos."
    )

    IJ.log("CSV por núcleo:")
    IJ.log(perNucleusOutput.getAbsolutePath())

    IJ.log("CSV por imagen:")
    IJ.log(perImageOutput.getAbsolutePath())

    IJ.log("================================================")

} catch (Exception exception) {

    IJ.log("ERROR guardando los resultados finales:")
    IJ.log(exception.toString())

    throw exception
}


// =====================================================
// FUNCIONES AUXILIARES
// =====================================================

/**
 * Convierte el nombre de una imagen en el nombre esperado
 * de su máscara de Cellpose.
 */
String createMaskFileName(String rawFileName) {

    String lowerName =
            rawFileName.toLowerCase(Locale.ROOT)

    if (lowerName.endsWith(".tiff")) {

        return rawFileName.substring(
                0,
                rawFileName.length() - 5
        ) + "_cp_masks.tif"
    }

    if (lowerName.endsWith(".tif")) {

        return rawFileName.substring(
                0,
                rawFileName.length() - 4
        ) + "_cp_masks.tif"
    }

    return rawFileName +
            "_cp_masks.tif"
}


/**
 * Elimina del citoplasma virtual:
 *
 * 1. Vóxeles fuera de la imagen.
 * 2. Vóxeles pertenecientes a cualquier núcleo.
 */
Object3DVoxels cleanVirtualCytoObject(
        Object3DVoxels inputObject,
        ImageInt nucleiLabels,
        int imageWidth,
        int imageHeight,
        int imageDepth,
        boolean excludeAllNuclei,
        double resolutionXY,
        double resolutionZ,
        String unit
) {

    if (inputObject == null) {
        return null
    }

    LinkedList<Voxel3D> validVoxels =
            new LinkedList<Voxel3D>()

    List<Voxel3D> inputVoxels =
            inputObject.getVoxels()

    if (inputVoxels == null ||
            inputVoxels.isEmpty()) {

        return null
    }

    for (Voxel3D voxel : inputVoxels) {

        int x = voxel.getRoundX()
        int y = voxel.getRoundY()
        int z = voxel.getRoundZ()

        /*
         * No incluir vóxeles fuera de la imagen.
         */
        if (x < 0 || x >= imageWidth ||
                y < 0 || y >= imageHeight ||
                z < 0 || z >= imageDepth) {

            continue
        }

        /*
         * No incluir ningún vóxel nuclear.
         */
        if (excludeAllNuclei) {

            int nucleusLabel =
                    nucleiLabels.getPixelInt(
                            x,
                            y,
                            z
                    )

            if (nucleusLabel > 0) {
                continue
            }
        }

        validVoxels.add(
                new Voxel3D(
                        x,
                        y,
                        z,
                        1
                )
        )
    }

    if (validVoxels.isEmpty()) {
        return null
    }

    Object3DVoxels cleanedObject =
            new Object3DVoxels(validVoxels)

    cleanedObject.setValue(1)

    cleanedObject.setCalibration(
            resolutionXY,
            resolutionZ,
            unit
    )

    return cleanedObject
}


/**
 * Comprueba las dimensiones XYZ.
 */
boolean sameDimensions(
        ImagePlus firstImage,
        ImagePlus secondImage
) {

    if (firstImage == null ||
            secondImage == null) {

        return false
    }

    return firstImage.getWidth() ==
            secondImage.getWidth() &&

            firstImage.getHeight() ==
            secondImage.getHeight() &&

            firstImage.getNSlices() ==
            secondImage.getNSlices()
}


/**
 * Devuelve las dimensiones de una imagen como texto.
 */
String dimensionsAsString(
        ImagePlus image
) {

    if (image == null) {
        return "null"
    }

    return image.getWidth() +
            " x " +
            image.getHeight() +
            " x " +
            image.getNSlices()
}


/**
 * Convierte Number/Object a double.
 */
double numberValue(Object value) {

    if (value == null) {
        return 0.0d
    }

    if (value instanceof Number) {

        return ((Number) value)
                .doubleValue()
    }

    return Double.parseDouble(
            value.toString()
    )
}


/**
 * Calcula la media de una lista.
 */
double meanValues(
        List<Double> values
) {

    if (values == null ||
            values.isEmpty()) {

        return 0.0d
    }

    double total = 0.0d
    int validCount = 0

    for (Double value : values) {

        if (value != null &&
                !Double.isNaN(value) &&
                !Double.isInfinite(value)) {

            total += value.doubleValue()
            validCount++
        }
    }

    if (validCount == 0) {
        return 0.0d
    }

    return total /
            ((double) validCount)
}


/**
 * Calcula la suma de una lista.
 */
double sumValues(
        List<Double> values
) {

    if (values == null ||
            values.isEmpty()) {

        return 0.0d
    }

    double total = 0.0d

    for (Double value : values) {

        if (value != null &&
                !Double.isNaN(value) &&
                !Double.isInfinite(value)) {

            total += value.doubleValue()
        }
    }

    return total
}