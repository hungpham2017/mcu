months = ['January', 'Febraury', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December']; 
var theDate = new Date(document.lastModified); 
theDate.setTime((theDate.getTime()) ) 
with (theDate) { 
document.write("<i>Copyright &copy; 2019, Hung Q. Pham. Last updated: " + months[getMonth()]+' '+ getDate() +', '+ getFullYear()+' '+ getHours()+':'+ getMinutes()+"</i>") 