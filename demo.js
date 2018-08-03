

const numeric = require('numeric');
const math=require('mathjs');

function RMatrixSolveLS(A,NRows,NCols,B,obj){
    var U=numeric.svd(A).U;
    var V=numeric.svd(A).V;
    var S=numeric.svd(A).S;
  
    //tạo ma trận Z là nghịch đảo của S
    var Z=new Array(S.length);
    for(var k=0;k<S.length;k++)
    {
        Z[k]=new Array(S.length);
    }
    for(var i=0;i<S.length;i++)
        for(var j=0;j<S.length;j++)
        {
            if(i===j)
            {
                Z[i][j]=1.0/S[i];
            }
            else
                Z[i][j]=0;
        }
    var UT=math.transpose(U);
    var VT=math.transpose(V);
    var result;
    if(NRows>=NCols)
    {
        result=math.multiply(V,Z);
        result=math.multiply(result,UT);
        result=math.multiply(result,B);
    }
    else
    {
        result=math.multiply(U,Z);
        result=math.multiply(result,VT);
        result=math.multiply(result,B);
    }

    return result;

}




var readline = require('readline-sync');
var input=new Array(16);
var name=new Array(10);
name[0]="ammonium orthomolybdate";
name[1]="Borid Acid";
name[2]="calcium nitrate tetrahydrate";
name[3]="copper nitrate hexahydrate";
name[4]="iron dtpa";
name[5]="magnesium nitrate hexahydrate";
name[6]="manganese sulfate monohydrate";
name[7]="phosphoric acid 30%";
name[8]="potassium nitrate";
name[9]="zinc sulfate dihydrate";
console.log("----------------------------------HYDRO BUDDY---------------------------------");
console.log("                    *********Console Version 1.0*********                       \n ");
console.log("                    * Thử nghiệm trên tập hợp các loại phân:\n                     + ammonium orthomolybdate\n                     + Borid Acid\n                     + calcium nitrate tetrahydrate\n                     + copper nitrate hexahydrate\n                     + iron dtpa\n                     + magnesium nitrate hexahydrate\n                     + manganese sulfate monohydrate\n                     + phosphoric acid 30%\n                     + potassium nitrate\n                     + zinc sulfate dihydrate\n")
console.log("                    * Xin mời nhập các giá trị dưỡng chất sau:\n");
var NO3 = readline.question("*********************Nhập hàm lượng N(NO3-):********************* \n");
input[0]=parseFloat(NO3);
var K = readline.question("*********************Nhập hàm lượng K:********************* \n");
input[1]=parseFloat(K);

var P = readline.question("*********************Nhập hàm lượng P:********************* \n");
input[2]=parseFloat(P);

var Mg = readline.question("*********************Nhập hàm lượng Mg:********************* \n");
input[3]=parseFloat(Mg);
var Ca = readline.question("*********************Nhập hàm lượng Ca:********************* \n");
input[4]=parseFloat(Ca);
var S = readline.question("*********************Nhập hàm lượng S:********************* \n");
input[5]=parseFloat(S);
var Fe = readline.question("*********************Nhập hàm lượng Fe:********************* \n");
input[6]=parseFloat(Fe);
var Zn = readline.question("*********************Nhập hàm lượng Zn:********************* \n");
input[7]=parseFloat(Zn);
var Bo = readline.question("*********************Nhập hàm lượng B:********************* \n");
input[8]=parseFloat(Bo);
var Cu = readline.question("*********************Nhập hàm lượng Cu:********************* \n");
input[9]=parseFloat(Cu);
var Mo = readline.question("*********************Nhập hàm lượng Mo:********************* \n");
input[10]=parseFloat(Mo);
var Na = readline.question("*********************Nhập hàm lượng Na:********************* \n");
input[11]=parseFloat(Na);
var Si = readline.question("*********************Nhập hàm lượng Si:********************* \n");
input[12]=parseFloat(Si);
var Cl = readline.question("*********************Nhập hàm lượng Cl:********************* \n");
input[13]=parseFloat(Cl);
var Mn = readline.question("*********************Nhập hàm lượng Mn:********************* \n");
input[14]=parseFloat(Mn);
var NH4 = readline.question("*********************Nhập hàm lượng N(NH4+):********************* \n");
input[15]=parseFloat(NH4);
var ec_contribution=new Array(16);
ec_contribution[0] = 71.46 / 14;
ec_contribution[1] = 106 / 24.30;
ec_contribution[2] = 57 / 31;
ec_contribution[3] = 119 / 40;
ec_contribution[4] = 160 / 32;
ec_contribution[5] = 108.0 / 56;
ec_contribution[6] = 0;
ec_contribution[7] = 0;
ec_contribution[8] = 0;
ec_contribution[9]= 50.01 / 23;
ec_contribution[10]= 100 / 28.09;
ec_contribution[11]= 76.35 / 35.5;
ec_contribution[12]= 0;
ec_contribution[13]= 73.5 / 14;
ec_contribution[14]= 0;
ec_contribution[15] = 73 / 39;

//tạo ma trận tính toán
var A=[[0,0,1.186,0.948,0,1.093,0,0,1.3856,0],
        [0,0,0,0,0,0,0,0,3.867,0],
        [0,0,0,0,0,0,0,0.9483,0,0],
        [0,0,0,0,0,0.948,0,0,0,0],
        [0,0,1.697,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,1.8973,0,0,1.624],
        [0,0,0,0,0.7,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,3.311],
        [0,1.7482,0,0,0,0,0,0,0,0],
        [4.8943,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,3.2504,0,0,0],
        [1.4294,0,0,0,0,0,0,0,0,0]];
        var A1=[[0,0,1.186,0.948,0,1.093,0,0,1.3856,0],
        [0,0,0,0,0,0,0,0,3.867,0],
        [0,0,0,0,0,0,0,0.9483,0,0],
        [0,0,0,0,0,0.948,0,0,0,0],
        [0,0,1.697,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,1.8973,0,0,1.624],
        [0,0,0,0,0.7,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,3.311],
        [0,1.7482,0,0,0,0,0,0,0,0],
        [0,0,0,2.149,0,0,0,0,0,0],
        [4.8943,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,3.2504,0,0,0],
        [1.4294,0,0,0,0,0,0,0,0,0]];
var B=new Array(15);


for(var i=0;i<15;i++)
{
    if(i<9)
    {
        B[i]=input[i];
    }
    else 
    {
        B[i]=input[i+1];
    }
}


var x=new Array(10);

x=RMatrixSolveLS(A,15,10,B);
console.log("\n                     ************KHỐI LƯỢNG TÍNH TOÁN***********");

     for(var i=0;i<10;i++)
     {
         if(x[i]<0)
         {
             console.log("                     **WARNING: Phân bón "+name[i]+" được giải ra giá trị âm, tạm thời sẽ được gán bằng 0, sau khi tính toán xong hãy điều chỉnh lại các loại phân bạn chọn.\n")
             x[i]=0;
         }
         x[i]=math.round(x[i]*1000)/1000;
         console.log("                     +   "+name[i]);
         console.log("                     +   "+x[i]+"gram\n");
     }
console.log("\n                     ****************TỔNG GIÁ TIỀN**************");
var money=0;
for(var i=0;i<10;i++)
{
    money+=x[i]*100/1000;
}
console.log("                     "+math.round(money*10)/10 + " USD\n\n")
var output=new Array(16);
for(var i=0;i<16;i++)
{
    output[i]=0;
    for(var j=0;j<10;j++)
    {
        output[i]+=A1[i][j]*x[j];
    }
    output[i]=math.round(output[i]*1000)/1000;
}
var predicted_ec=0;
for(var i=0;i<16;i++)
{
    predicted_ec+=output[i]*ec_contribution[i];
}
predicted_ec=math.round((0.65*predicted_ec/1000)*10)/10;
console.log("\n                     ****************HÀM LƯỢNG EC DỰ KIẾN**************");
console.log("                     *  EC= "+ predicted_ec +" mS/cm \n\n")
var error=new Array(16);
for(var i=0;i<16;i++)
{
    if(input[i]==0)
    {
        error[i]=0;
    }
    else
    {
        error[i]=((output[i]-input[i])/input[i])*100;
        error[i]=math.round(error[i]*10)/10;
    }
}


console.log("**********************************KẾT QUẢ*********************************\n")
console.log("**************************************************************************");
console.log("*     CHẤT     *       LT(ppm)      *   KQ(ppm)  *      error(%)   *");
console.log("**************************************************************************");
console.log("*   N(NO3-)    *       "+input[0]+"       *     "+output[0]+"       *     "+error[0]+"    *");
console.log("*      K       *       "+input[1]+"       *     "+output[1]+"       *     "+error[1]+"    *");
console.log("*      P       *       "+input[2]+"       *     "+output[2]+"       *     "+error[2]+"    *");
console.log("*      Mg      *       "+input[3]+"       *     "+output[3]+"       *     "+error[3]+"    *");
console.log("*      Ca      *       "+input[4]+"       *     "+output[4]+"       *     "+error[4]+"    *");
console.log("*      S       *       "+input[5]+"       *     "+output[5]+"       *     "+error[5]+"    *");
console.log("*      Fe      *       "+input[6]+"       *     "+output[6]+"       *     "+error[6]+"    *");
console.log("*      Zn      *       "+input[7]+"       *     "+output[7]+"       *     "+error[7]+"    *");
console.log("*      B       *       "+input[8]+"       *     "+output[8]+"       *     "+error[8]+"    *");
console.log("*      Cu      *       "+input[9]+"       *     "+output[9]+"       *     "+error[9]+"    *");
console.log("*      Mo      *       "+input[10]+"       *     "+output[10]+"       *     "+error[10]+"    *");
console.log("*      Na      *       "+input[11]+"       *     "+output[11]+"       *     "+error[11]+"    *");
console.log("*      Si      *       "+input[12]+"       *     "+output[12]+"       *     "+error[12]+"    *");
console.log("*      Cl      *       "+input[13]+"       *     "+output[13]+"       *     "+error[13]+"    *");
console.log("*      Mn      *       "+input[14]+"       *     "+output[14]+"       *     "+error[14]+"    *");
console.log("*   N(NH4+)    *       "+input[15]+"       *     "+output[15]+"       *     "+error[15]+"    *");
console.log("**************************************************************************");

console.log("\n\nDựa vào tên của các loại phân đã biết, ta có thể chia thành số bồn tối thiểu là 2 bồn \n");
console.log("******Trong đó bồn A gồm những chất sau******");
console.log(name[2]+"\n"+name[4]+"\n"+name[5]+"\n"+name[8]+"\n");
console.log("******Trong đó bồn B gồm những chất sau******");
console.log(name[0]+"\n"+name[1]+"\n"+name[3]+"\n"+name[6]+"\n"+name[7]+"\n"+name[9]+"\n");
var sobon_nhap = readline.question("************Còn bạn muốn sử dụng bao nhiêu bồn?*********\n");
var sobon=parseInt(sobon_nhap);
var option = readline.question("*********************Bạn muốn nhập theo option nào(1/2) ?\n    1.Theo số cây\n    2.Theo mật độ cây trên hecta và diện tích trồng\n ");
var opt=parseInt(option);
while(opt!=1&&opt!=2)
{
    console.log("Cú pháp đã sai. Mời bạn nhập lại.\n")
    var option_temp = readline.question("*********************Bạn muốn nhập theo option nào(1/2) ?\n    1.Theo số cây\n    2.Theo mật độ cây trên hecta và diện tích trồng\n ");
    opt=parseInt(option_temp);
}
if(opt==1)
{
    var numTree = readline.question("*********************Nhập số cây:********************* \n");
    var soCay=parseInt(numTree);
}
else if(opt==2)
{
    var numTree_ha = readline.question("*********************Mật độ cây trên hecta (số cây/hecta):********************* \n");
    var soCay_ha=parseFloat(numTree_ha);
    var dientich = readline.question("*********************Nhập diện tích trồng:********************* \n");
    var dienTichTrong=parseFloat(dientich);
    var soCay1=parseInt(dienTichTrong/soCay_ha);
}


var percent=new Array(10);
var phan0 = readline.question("*********************Nhập % so với độ bão hòa của phân "+name[0]+"\n");
percent[0]=parseFloat(phan0);
var phan1 = readline.question("*********************Nhập % so với độ bão hòa của phân "+name[1]+"\n");
percent[1]=parseFloat(phan1);
var phan2 = readline.question("*********************Nhập % so với độ bão hòa của phân "+name[2]+" \n");
percent[2]=parseFloat(phan2);
var phan3 = readline.question("*********************Nhập % so với độ bão hòa của phân "+name[3]+" \n");
percent[3]=parseFloat(phan3);
var phan4 = readline.question("*********************Nhập % so với độ bão hòa của phân "+name[4]+" \n");
percent[4]=parseFloat(phan4);
var phan5 = readline.question("*********************Nhập % so với độ bão hòa của phân "+name[5]+" \n");
percent[5]=parseFloat(phan5);
var phan6 = readline.question("*********************Nhập % so với độ bão hòa của phân "+name[6]+" \n");
percent[6]=parseFloat(phan6);
var phan7 = readline.question("*********************Nhập % so với độ bão hòa của phân "+name[7]+"\n");
percent[7]=parseFloat(phan7);
var phan8 = readline.question("*********************Nhập % so với độ bão hòa của phân "+name[8]+" \n");
percent[8]=parseFloat(phan8);
var phan9 = readline.question("*********************Nhập % so với độ bão hòa của phân "+name[9]+"\n");
percent[9]=parseFloat(phan9);
var litMotBon=100.0;
var baohoa=new Array(10);//nồng độ bão hòa gốc
baohoa[0]=200.1;
baohoa[1]=200.1;
baohoa[2]=200.1;
baohoa[3]=200.1;
baohoa[4]=200.1;
baohoa[5]=200.1;
baohoa[6]=200.1;
baohoa[7]=200.1;
baohoa[8]=200.1;
baohoa[9]=200.1;
var nongdo=new Array(10) //Nồng độ thực tế
for(var i=0;i<10;i++)
{
    nongdo[i]=percent[i]*baohoa[i]/100;
}
//tính số lượng phân cần
var khoiLuongCan=new Array(10);
if(opt==1)//dùng socay
{
    for(var i=0;i<10;i++)
    {
        khoiLuongCan[i]=x[i]*soCay;
    }
}
else//dùng socay1
{
    for(var i=0;i<10;i++)
    {
        khoiLuongCan[i]=x[i]*soCay1;
    }
}
var luongNuocCan=new Array(10);//tính được lượng nước cần cho mỗi loại cây
for (var i=0;i<10;i++)
{
    luongNuocCan[i]=((khoiLuongCan[i]*100)/nongdo[i])/1000;
    luongNuocCan[i]=math.round(luongNuocCan[i]*1000)/1000;
}



//xuất trường hợp 2 bồn
if(sobon==2)
{
    console.log("***********Bạn đã chọn pha 2 bồn***********\n\n")
    var bonA=new Array();
    bonA.push(luongNuocCan[2]);
    bonA.push(luongNuocCan[4]);
    bonA.push(luongNuocCan[5]);
    bonA.push(luongNuocCan[8]);
    var bonB=new Array();
    bonB.push(luongNuocCan[0]);
    bonB.push(luongNuocCan[1]);
    bonB.push(luongNuocCan[3]);
    bonB.push(luongNuocCan[6]);
    bonB.push(luongNuocCan[7]);
    bonB.push(luongNuocCan[9]);
    var nuocA=math.max(bonA);
    var nuocB=math.max(bonB);
    console.log("***********Phân ở bồn A là:\n "+name[2]+" : "+khoiLuongCan[2]/1000+"kg\n "+name[4]+" : "+khoiLuongCan[4]/1000+"kg\n"+name[5]+" : "+khoiLuongCan[5]/1000+"kg\n"+name[8]+" : "+khoiLuongCan[8]/1000+"kg\n Với lượng nước ở bồn A là: "+nuocA+" lít ");
    console.log("***********Phân ở bồn B là:\n "+name[0]+" : "+khoiLuongCan[0]/1000+"kg\n "+name[1]+" : "+khoiLuongCan[1]/1000+"kg\n"+name[3]+" : "+khoiLuongCan[3]/1000+"kg\n"+name[6]+" : "+khoiLuongCan[6]/1000+"kg\n"+name[7]+" : "+khoiLuongCan[7]/1000+"kg\n"+name[9]+" : "+khoiLuongCan[9]/1000+"kg\n Với lượng nước ở bồn B là: "+nuocB+" lít ");

}
else if(sobon==3)
{
    console.log("***********Bạn đã chọn pha 3 bồn***********\n\n")
    var bonA=new Array();
    bonA.push(luongNuocCan[2]);
    bonA.push(luongNuocCan[4]);
    bonA.push(luongNuocCan[5]);
    bonA.push(luongNuocCan[8]);
    var bonB=new Array();
    bonB.push(luongNuocCan[0]);
    bonB.push(luongNuocCan[1]);
    bonB.push(luongNuocCan[3]);
    bonB.push(luongNuocCan[6]);
    bonB.push(luongNuocCan[7]);
    bonB.push(luongNuocCan[9]);
    var bonC=new Array();
    var indexMaxA1 = bonA.indexOf(math.max(bonA));//index của giá trị lớn nhất của A
    //tìm index của đứa lớn nhất trong luongnuoccan để loại nó đi khi xuất ra
    var indexmain=luongNuocCan.indexOf(math.max(bonA));
    bonC.push(bonA[indexMaxA1]);//bồn C giờ chứa chất có lượng nước max từ A
    if (indexMaxA1 > -1) {
        bonA.splice(indexMaxA1, 1);
    }
    var nuocA=math.max(bonA);
    var nuocB=math.max(bonB);
    var nuocC=bonC[0];
    console.log("Phân ở bồn A là:");
    for(var i=0;i<10;i++)
    {
        if((i==2||i==4||i==5||i==8)&&(i!=indexmain))
        {
            console.log(name[i]+" : "+khoiLuongCan[i]/1000+" kg");
        }
    }
    console.log("Với lượng nước bồn A là: "+nuocA+" lít");

    console.log("Phân ở bồn B là:");
    for(var i=0;i<10;i++)
    {
        if(i==0||i==1||i==3||i==6||i==7||i==9)
        {
            console.log(name[i]+" : "+khoiLuongCan[i]/1000+" kg");
        }
    }
    console.log("Với lượng nước bồn B là: "+nuocB+" lít");


    console.log("Phân ở bồn C là:");
    console.log(name[indexmain]+" : "+khoiLuongCan[indexmain]/1000+" kg");
    console.log("Với lượng nước bồn C là: "+nuocC+" lít");

}
else if(sobon==4)
{
    console.log("***********Bạn đã chọn pha 4 bồn***********\n\n")
    var bonA=new Array();
    bonA.push(luongNuocCan[2]);
    bonA.push(luongNuocCan[4]);
    bonA.push(luongNuocCan[5]);
    bonA.push(luongNuocCan[8]);
    var bonB=new Array();
    bonB.push(luongNuocCan[0]);
    bonB.push(luongNuocCan[1]);
    bonB.push(luongNuocCan[3]);
    bonB.push(luongNuocCan[6]);
    bonB.push(luongNuocCan[7]);
    bonB.push(luongNuocCan[9]);
    var bonC=new Array();
    var bonD=new Array();
    var indexMaxA1 = bonA.indexOf(math.max(bonA));//index của giá trị lớn nhất của A
    var indexMaxB1 = bonB.indexOf(math.max(bonB));//index của giá trị lớn nhất của B
    var indexmain1=luongNuocCan.indexOf(math.max(bonA));
    var indexmain2=luongNuocCan.indexOf(math.max(bonB));
    bonC.push(bonA[indexMaxA1]);//bồn C giờ chứa chất có lượng nước max từ A
    if (indexMaxA1 > -1) {
        bonA.splice(indexMaxA1, 1);
    }
    bonD.push(bonB[indexMaxB1]);//bồn D giờ chứa chất có lượng nước max từ B
    if (indexMaxB1 > -1) {
        bonB.splice(indexMaxB1, 1);
    }

    var nuocA=math.max(bonA);
    var nuocB=math.max(bonB);
    var nuocC=bonC[0];
    var nuocD=bonD[0];
    console.log("Phân ở bồn A là:");
    for(var i=0;i<10;i++)
    {
        if((i==2||i==4||i==5||i==8)&&(i!=indexmain1))
        {
            console.log(name[i]+" : "+khoiLuongCan[i]/1000+" kg");
        }
    }
    console.log("Với lượng nước bồn A là: "+nuocA+" lít");

    console.log("Phân ở bồn B là:");
    for(var i=0;i<10;i++)
    {
        if((i==0||i==1||i==3||i==6||i==7||i==9)&&(i!=indexmain2))
        {
            console.log(name[i]+" : "+khoiLuongCan[i]/1000+" kg");
        }
    }
    console.log("Với lượng nước bồn B là: "+nuocB+" lít");


    console.log("Phân ở bồn C là:");
    console.log(name[indexmain1]+" : "+khoiLuongCan[indexmain1]/1000+" kg");
    console.log("Với lượng nước bồn C là: "+nuocC+" lít");
    console.log("Phân ở bồn D là:");
    console.log(name[indexmain2]+" : "+khoiLuongCan[indexmain2]/1000+" kg");
    console.log("Với lượng nước bồn D là: "+nuocD+" lít");

}
else if(sobon==5)
{
    console.log("***********Bạn đã chọn pha 5 bồn***********\n\n")
    var bonA=new Array();
    bonA.push(luongNuocCan[2]);
    bonA.push(luongNuocCan[4]);
    bonA.push(luongNuocCan[5]);
    bonA.push(luongNuocCan[8]);
    var bonB=new Array();
    bonB.push(luongNuocCan[0]);
    bonB.push(luongNuocCan[1]);
    bonB.push(luongNuocCan[3]);
    bonB.push(luongNuocCan[6]);
    bonB.push(luongNuocCan[7]);
    bonB.push(luongNuocCan[9]);
    var bonC=new Array();
    var bonD=new Array();
    var bonE=new Array();
    var indexMaxA1 = bonA.indexOf(math.max(bonA));//index của giá trị lớn nhất của A
    var indexMaxB1 = bonB.indexOf(math.max(bonB));//index của giá trị lớn nhất của B
    var indexmain1=luongNuocCan.indexOf(math.max(bonA));
    var indexmain2=luongNuocCan.indexOf(math.max(bonB));
    bonC.push(bonA[indexMaxA1]);//bồn C giờ chứa chất có lượng nước max từ A
    if (indexMaxA1 > -1) {
        bonA.splice(indexMaxA1, 1);
    }
    var indexmain3=luongNuocCan.indexOf(math.max(bonA));
    bonD.push(bonB[indexMaxB1]);//bồn D giờ chứa chất có lượng nước max từ B
    if (indexMaxB1 > -1) {
        bonB.splice(indexMaxB1, 1);
    }
    var indexMaxA2=bonA.indexOf(math.max(bonA));
    bonE.push(bonA[indexMaxA2]);
    if (indexMaxA2 > -1) {//xóa thằng cao nhì
        bonA.splice(indexMaxA2, 1);
    }
    var nuocA=math.max(bonA);
    var nuocB=math.max(bonB);
    var nuocC=bonC[0];
    var nuocD=bonD[0];
    var nuocE=bonE[0];
    console.log("Phân ở bồn A là:");
    for(var i=0;i<10;i++)
    {
        if((i==2||i==4||i==5||i==8)&&(i!=indexmain1)&&(i!=indexmain3))
        {
            console.log(name[i]+" : "+khoiLuongCan[i]/1000+" kg");
        }
    }
    console.log("Với lượng nước bồn A là: "+nuocA+" lít\n");

    console.log("Phân ở bồn B là:");
    for(var i=0;i<10;i++)
    {
        if((i==0||i==1||i==3||i==6||i==7||i==9)&&(i!=indexmain2))
        {
            console.log(name[i]+" : "+khoiLuongCan[i]/1000+" kg");
        }
    }
    console.log("Với lượng nước bồn B là: "+nuocB+" lít\n");


    console.log("Phân ở bồn C là:");
    console.log(name[indexmain1]+" : "+khoiLuongCan[indexmain1]/1000+" kg");
    console.log("Với lượng nước bồn C là: "+nuocC+" lít\n");
    console.log("Phân ở bồn D là:");
    console.log(name[indexmain2]+" : "+khoiLuongCan[indexmain2]/1000+" kg");
    console.log("Với lượng nước bồn D là: "+nuocD+" lít\n");
    console.log("Phân ở bồn E là:");
    console.log(name[indexmain3]+" : "+khoiLuongCan[indexmain3]/1000+" kg");
    console.log("Với lượng nước bồn E là: "+nuocE+" lít");    
}
else if(sobon==6)
{
    console.log("***********Bạn đã chọn pha 5 bồn***********\n\n")
    var bonA=new Array();
    bonA.push(luongNuocCan[2]);
    bonA.push(luongNuocCan[4]);
    bonA.push(luongNuocCan[5]);
    bonA.push(luongNuocCan[8]);
    var bonB=new Array();
    bonB.push(luongNuocCan[0]);
    bonB.push(luongNuocCan[1]);
    bonB.push(luongNuocCan[3]);
    bonB.push(luongNuocCan[6]);
    bonB.push(luongNuocCan[7]);
    bonB.push(luongNuocCan[9]);
    var bonC=new Array();
    var bonD=new Array();
    var bonE=new Array();
    var bonF=new Array();
    var indexMaxA1 = bonA.indexOf(math.max(bonA));//index của giá trị lớn nhất của A
    var indexMaxB1 = bonB.indexOf(math.max(bonB));//index của giá trị lớn nhất của B
    var indexmain1=luongNuocCan.indexOf(math.max(bonA));
    var indexmain2=luongNuocCan.indexOf(math.max(bonB));
    bonC.push(bonA[indexMaxA1]);//bồn C giờ chứa chất có lượng nước max từ A
    if (indexMaxA1 > -1) {
        bonA.splice(indexMaxA1, 1);
    }
    var indexmain3=luongNuocCan.indexOf(math.max(bonA));
    bonD.push(bonB[indexMaxB1]);//bồn D giờ chứa chất có lượng nước max từ B
    if (indexMaxB1 > -1) {
        bonB.splice(indexMaxB1, 1);
    }
    var indexmain4=luongNuocCan.indexOf(math.max(bonB));
    var indexMaxA2=bonA.indexOf(math.max(bonA));
    var indexMaxB2=bonB.indexOf(math.max(bonB));
    bonE.push(bonA[indexMaxA2]);
    if (indexMaxA2 > -1) {//xóa thằng cao nhì
        bonA.splice(indexMaxA2, 1);
    }
    bonF.push(bonB[indexMaxB2]);
    if (indexMaxB2 > -1) {//xóa thằng cao nhì
        bonB.splice(indexMaxB2, 1);
    }
    var nuocA=math.max(bonA);
    var nuocB=math.max(bonB);
    var nuocC=bonC[0];
    var nuocD=bonD[0];
    var nuocE=bonE[0];
    var nuocF=bonF[0];
    console.log("Phân ở bồn A là:");
    for(var i=0;i<10;i++)
    {
        if((i==2||i==4||i==5||i==8)&&(i!=indexmain1)&&(i!=indexmain3))
        {
            console.log(name[i]+" : "+khoiLuongCan[i]/1000+" kg");
        }
    }
    console.log("Với lượng nước bồn A là: "+nuocA+" lít\n");

    console.log("Phân ở bồn B là:");
    for(var i=0;i<10;i++)
    {
        if((i==0||i==1||i==3||i==6||i==7||i==9)&&(i!=indexmain2)&&(i!=indexmain4))
        {
            console.log(name[i]+" : "+khoiLuongCan[i]/1000+" kg");
        }
    }
    console.log("Với lượng nước bồn B là: "+nuocB+" lít\n");
    console.log("Phân ở bồn C là:");
    console.log(name[indexmain1]+" : "+khoiLuongCan[indexmain1]/1000+" kg");
    console.log("Với lượng nước bồn C là: "+nuocC+" lít\n");
    console.log("Phân ở bồn D là:");
    console.log(name[indexmain2]+" : "+khoiLuongCan[indexmain2]/1000+" kg");
    console.log("Với lượng nước bồn D là: "+nuocD+" lít\n");
    console.log("Phân ở bồn E là:");
    console.log(name[indexmain3]+" : "+khoiLuongCan[indexmain3]/1000+" kg");
    console.log("Với lượng nước bồn E là: "+nuocE+" lít\n");
    console.log("Phân ở bồn F là:");
    console.log(name[indexmain4]+" : "+khoiLuongCan[indexmain4]/1000+" kg");
    console.log("Với lượng nước bồn F là: "+nuocF+" lít\n");    
}
//3 bồn: bỏ chất lít cao nhất của bồn A qua bồn C
//4 bồn: bỏ chất lít cao nhất của A qua C, B qua D
//5 bồn: bỏ chất lít cao nhất của A qua C, B cao nhất qua D, cao nhì qua E