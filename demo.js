

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
console.log("-------------------------HYDRO BUDDY------------------------");
console.log("                    ***Console Version 1.0***                       \n ");
console.log("   * Thử nghiệm trên tập hợp các loại phân:\n    + ammonium orthomolybdate\n    + Borid Acid\n    + calcium nitrate tetrahydrate\n    + copper nitrate hexahydrate\n    + iron dtpa\n    + magnesium nitrate hexahydrate\n    + manganese sulfate monohydrate\n    + phosphoric acid 30%\n    + potassium nitrate\n    + zinc sulfate dihydrate\n")
console.log("   * Xin mời nhập các giá trị dưỡng chất sau:\n");
var NO3 = readline.question("***Nhập hàm lượng N(NO3-):*** \n");
input[0]=parseFloat(NO3);
var NH4 = readline.question("***Nhập hàm lượng N(NH4+):*** \n");
input[1]=parseFloat(NH4);
var P = readline.question("***Nhập hàm lượng P:*** \n");
input[2]=parseFloat(P);
var K = readline.question("***Nhập hàm lượng K:*** \n");
input[3]=parseFloat(K);
var Mg = readline.question("***Nhập hàm lượng Mg:*** \n");
input[4]=parseFloat(Mg);
var Ca = readline.question("***Nhập hàm lượng Ca:*** \n");
input[5]=parseFloat(Ca);
var S = readline.question("***Nhập hàm lượng S:*** \n");
input[6]=parseFloat(S);
var Fe = readline.question("***Nhập hàm lượng Fe:*** \n");
input[7]=parseFloat(Fe);
var Zn = readline.question("***Nhập hàm lượng Zn:*** \n");
input[8]=parseFloat(Zn);
var Bo = readline.question("***Nhập hàm lượng B:*** \n");
input[9]=parseFloat(Bo);
var Mn = readline.question("***Nhập hàm lượng Mn:*** \n");
input[10]=parseFloat(Mn);
var Cu = readline.question("***Nhập hàm lượng Cu:*** \n");
input[11]=parseFloat(Cu);
var Mo = readline.question("***Nhập hàm lượng Mo:*** \n");
input[12]=parseFloat(Mo);
var Na = readline.question("***Nhập hàm lượng Na:*** \n");
input[13]=parseFloat(Na);
var Si = readline.question("***Nhập hàm lượng Si:*** \n");
input[14]=parseFloat(Si);
var Cl = readline.question("***Nhập hàm lượng Cl:*** \n");
input[15]=parseFloat(Cl);


//tạo ma trận tính toán
var A=[[0,0,1.186,0.948,0,1.093,0,0,1.3856,0],
        [1.4294,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0.9483,0,0],
        [0,0,0,0,0,0,0,0,3.867,0],
        [0,0,0,0,0,0.948,0,0,0,0],
        [0,0,1.697,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,1.8973,0,0,1.624],
        [0,0,0,0,0.7,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,3.311],
        [0,1.7482,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,3.2504,0,0,0],
        [0,0,0,2.149,0,0,0,0,0,0],
        [4.8943,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0]];
var B=new Array(16);


for(var i=0;i<16;i++)
{
    B[i]=input[i];
}

console.log("\n************KHỐI LƯỢNG TÍNH TOÁN***********");
var x=new Array(10);
x=RMatrixSolveLS(A,16,10,B);
     for(var i=0;i<10;i++)
     {
         if(x[i]<0)
         {
             x[i]=0;
         }
         x[i]=math.round(x[i]*1000)/1000;
         console.log(name[i]);
         console.log(x[i]+"gram\n");
     }
console.log("\n****************TỔNG GIÁ TIỀN**************");
var money=0;
for(var i=0;i<10;i++)
{
    money+=x[i]*100/1000;
}
console.log(math.round(money*10)/10 + " USD\n\n")
var output=new Array(16);
for(var i=0;i<16;i++)
{
    output[i]=0;
    for(var j=0;j<10;j++)
    {
        output[i]+=A[i][j]*x[j];
    }
    output[i]=math.round(output[i]*1000)/1000;
}

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


console.log("**********************KẾT QUẢ***************************\n")
console.log("********************************************************");
console.log("*   CHẤT   *   LT(ppm)  *   KQ(ppm)  *      error(%)   *");
console.log("********************************************************");
console.log("* N(NO3-)  *   "+input[0]+"   *   "+output[0]+"     *     "+error[0]+"    *");
console.log("* N(NH4+)  *   "+input[1]+"   *   "+output[1]+"     *     "+error[1]+"    *");
console.log("*    P     *   "+input[2]+"   *   "+output[2]+"     *     "+error[2]+"    *");
console.log("*    K     *   "+input[3]+"   *   "+output[3]+"     *     "+error[3]+"    *");
console.log("*    Mg    *   "+input[4]+"   *   "+output[4]+"     *     "+error[4]+"    *");
console.log("*    Ca    *   "+input[5]+"   *   "+output[5]+"     *     "+error[5]+"    *");
console.log("*    S     *   "+input[6]+"   *   "+output[6]+"     *     "+error[6]+"    *");
console.log("*    Fe    *   "+input[7]+"   *   "+output[7]+"     *     "+error[7]+"    *");
console.log("*    Zn    *   "+input[8]+"   *   "+output[8]+"     *     "+error[8]+"    *");
console.log("*    B     *   "+input[9]+"   *   "+output[9]+"     *     "+error[9]+"    *");
console.log("*    Mn    *   "+input[10]+"   *   "+output[10]+"     *     "+error[10]+"    *");
console.log("*    Cu    *   "+input[11]+"   *   "+output[11]+"     *     "+error[11]+"    *");
console.log("*    Mo    *   "+input[12]+"   *   "+output[12]+"     *     "+error[12]+"    *");
console.log("*    Na    *   "+input[13]+"   *   "+output[13]+"     *     "+error[13]+"    *");
console.log("*    Si    *   "+input[14]+"   *   "+output[14]+"     *     "+error[14]+"    *");
console.log("*    Cl    *   "+input[15]+"   *   "+output[15]+"     *     "+error[15]+"    *");
console.log("********************************************************");

