## 5.0.2 - 2026-08-04
- [GRD-1175](https://jira.oicr.on.ca/browse/GRD-1175) 
- Updated to support hg38_noAlt reference 
- Fixed an error in contig splitting 

## 5.0.1 - 2023-12-01
- [GRD-683](https://jira.oicr.on.ca/browse/GRD-683) 
- Added "readGroup" argument
- Added parameter "skipcheck", works on GBS-1964 branch years ago
- Removed using of picard tool https://gatk.broadinstitute.org/hc/en-us/articles/360037226472-AddOrReplaceReadGroups-Picard- from master branch, which is a different approach to address readGroup info issue.