class UTILS {
    // Function to remove Nextflow version from pipeline_software_mqc_versions.yml
    public static Object removeNextflowVersion(pipeline_software_mqc_versions) {
        def softwareVersions = path(pipeline_software_mqc_versions).yaml
        if (softwareVersions.containsKey("Workflow")) {
            softwareVersions.Workflow.remove("Nextflow")
        }
        return softwareVersions
    }

    // Recursively list all files in a directory and its sub-directories, matching a given suffix
    // TODO: use regex pattern instead of suffix?
    public static getAllFilesFromDir(dir, suffix) {
        def output = []
        new File(dir).eachFileRecurse() {
            if (it.name.toString().endsWith(suffix)) {
                output.add(it)
            }
        }
        return output.sort()
    }
    // Recursively list all files names in a directory and its sub-directories, matching a given suffix, return file names
    public static getAllFileNamesFromDir(dir, suffix) {
        def output = []
        new File(dir).eachFileRecurse() {
            if (it.name.toString().endsWith(suffix)) {
                output.add(it.toString().split("/")[-1])
            }
        }
        return output.sort()
    }

    // Recursively list all files names in a directory and its sub-directories, matching a given suffix, return if check if given string is in file
    public static checkAllFilesNamesFromDirForString(dir, suffix, string) {
        def output = []
        new File(dir).eachFileRecurse() {
            if (it.name.toString().endsWith(suffix)) {
                output.add(it.text.contains(string))
            }
        }
        return output.sort()
    }
}
