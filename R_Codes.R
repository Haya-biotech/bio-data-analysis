dir.create("clean_data")
dir.create("scripts")
dir.create("result")
cholesterol <- 230

if(cholesterol > 240){
  print("High Cholesterol")
}

Systolic_bp <- 130
 if(Systolic_bp< 120){
   print("Blood Pressure is normal")
 }else{
   print("Blood Pressure is high")
 }
####patient info ####
data<- read.csv(file.choose())
patient<- data
View(patient)
str(patient)
factor_cols<- c("gender","diagnosis", "smoker")

for( col in factor_cols){
  patient[[col]] <- as.factor(patient[[col]])
}
str((patient))

binary_cols<-('smoker')

for(cols in binary_cols){
  patient[[cols]]<- ifelse(patient$smoker =="Yes", 1, 0)
}
str(data)
str(patient)
write.csv(patient, file = "clean_data/patient_info.csv")

#### metadata file####
mdata<- read.csv(file.choose())
meta_data<-mdata
View(meta_data)
str(meta_data)
factor_cols2<-c("height", "gender")
for(coll in factor_cols2){
  meta_data[[coll]]<- as.factor(meta_data[[coll]])
}
str(meta_data)
binary_cols2<-("gender")

for(i in binary_cols2){
  meta_data[[i]]<- ifelse(meta_data$gender == "Female", 1, 0)
}
str(mdata)
str(meta_data)
write.csv(meta_data, file = "clean_data/meta_data.csv")

