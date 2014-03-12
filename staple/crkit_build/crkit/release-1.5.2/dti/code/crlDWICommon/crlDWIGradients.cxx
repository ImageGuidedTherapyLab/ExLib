#include "crlDWIGradients.h"
#include "crlDWIGradientDirectionsDef.h"
#include <time.h>

namespace crl {
namespace DWI {

/**********************************************************************************************//**
 * \fn	DWIGradients::DWIGradients( void )
 *
 * \brief	Default constructor. 
 *
 * \author	Benoit Scherrer
 * \date	August 2010
*************************************************************************************************/
DWIGradients::DWIGradients( void ):
m_NominalB0(0),
m_Verbose(false)
{
	m_DiffusionVectors = GradientDirectionContainerType::New();

}

/**********************************************************************************************//**
 * \fn	DWIGradients::~DWIGradients()
 *
 * \brief	Destructor. 
 *
 * \author	Benoit Scherrer
 * \date	August 2010
*************************************************************************************************/
DWIGradients::~DWIGradients()
{
}

/**********************************************************************************************//**
 * \fn	void DWIGradients::ClearAllGradientDirections()
 *
 * \brief	Clears all gradient directions. 
 *
 * \author	Benoit Scherrer
 * \date	August 2010
*************************************************************************************************/
void DWIGradients::ClearAllGradientDirections()
{
	m_DiffusionVectors->Initialize();

}

/**********************************************************************************************//**
 * \fn	DWIGradients::GradientSchemeITKType* DWIGradients::GetSchemeAsITKContainer() const
 *
 * \brief	Gets the gradient scheme as an itk container (GradientSchemeITKType). 
 *
 * \author	Benoit Scherrer
 * \date	August 2010
 *
 * \return	null if it fails, else the scheme as itk container. 
*************************************************************************************************/
DWIGradients::GradientSchemeITKType* DWIGradients::GetSchemeAsITKContainer() const
{ 
	return m_DiffusionVectors.GetPointer();
}

/**********************************************************************************************//**
 * \fn	DWIGradients::GradientSchemeStdType DWIGradients::GetSchemeAsStdVector() const
 *
 * \brief	Gets the scheme as an std::vector (GradientSchemeStdType). 
 *
 * \author	Benoit Scherrer
 * \date	August 2010
 *
 * \return	The scheme as std vector. 
*************************************************************************************************/
DWIGradients::GradientSchemeStdType DWIGradients::GetSchemeAsStdVector() const
{
	GradientSchemeStdType vector;

	for ( unsigned int i=0; i<m_DiffusionVectors->Size() ; i++ )
	{
		GradientDirectionType vect3d = m_DiffusionVectors->GetElement(i);
		vector.push_back( vect3d );
	}

	return vector;
}

void DWIGradients::ExtractNormalizedGradientsAndBValues(GradientSchemeStdType& grads, std::vector<double>& bvalues, bool includeB0 )
{
	if ( m_NominalB0==0 )
		throw itk::ExceptionObject(__FILE__,__LINE__,"INTERNAL ERROR. The nominal b-value must be set before calling ExtractNormalizedGradientsAndBValues", "");

	grads.clear();
	bvalues.clear();

	//-----------------------------------------------
	// For all gradients
	//-----------------------------------------------
	for ( unsigned int i=0; i<m_DiffusionVectors->Size() ; i++ )
	{
		GradientDirectionType vect3d = m_DiffusionVectors->GetElement(i);

		//-----------------------------------------------
		// Null gradient? 
		//-----------------------------------------------
		if ( vect3d.one_norm() <= 1e-8 )
		{
			if ( includeB0 )
			{
				vect3d[0]=vect3d[1]=vect3d[2]=0;
				grads.push_back(vect3d);
				bvalues.push_back(0);
			}
		}

		//-----------------------------------------------
		// Non null gradient
		//-----------------------------------------------
		else
		{
			double nn = vect3d.two_norm();
			double b = m_NominalB0 * nn * nn;

			grads.push_back( vect3d / nn );
			bvalues.push_back(b);
		}
	}
}


/**********************************************************************************************//**
 * \fn	void DWIGradients::GetStringDescriptionHelp(std::string& shortHelp, std::string& longHelp)
 *
 * \brief	Get the help (for TCLAP)
 *
 * \author	Benoit Scherrer
 * \date	August 2010
 *
 * \param [in,out]	shortHelp	The short help. 
 * \param [in,out]	longHelp	The long help. 
*************************************************************************************************/
void DWIGradients::GetStringDescriptionHelp(std::string& shortHelp, std::string& longHelp)
{
	shortHelp = "[NominalB0] ['0' Nex] ['S' NbDir BVal] ['C' 'Nex' '2/3'] ..." ;
	longHelp = "Describe the gradient scheme. \n The first element is always the nominal B0. \
			   Then several occurences of the commands '0', 'S' or 'C' can be defined: \n \
			   '0' [Nb B0] : number of B0 (by default add 1 if not defined) \n \
			   'S' [Number Directions] [b-value] : add a shell of N-directions at a b-value\n \
			   'C' [Number Repetition] [2/3] : add n-repetition of a cube corner of norm 2 or 3 " ;

}

/**********************************************************************************************//**
 * \fn	void DWIGradients::SetFromString(const std::string& gradientScheme)
 *
 * \brief	Sets the whole gradient scheme from a string description. 
 *
 * [NominalB0]  
 * 'S' [Number Directions] [b-value] : add a shell of N-directions at a b-value
 * 'C' [Number Repetition] [2/3] : add n-repetition of a cube corner of norm 2 or 3
 * '0' [Nb B0] : number of B0 (by default add 1).
 *
 * \author	Benoit Scherrer
 * \date	August 2010
 *
 * \param	gradientScheme	The gradient scheme. 
*************************************************************************************************/
void DWIGradients::SetFromString(const std::string& gradientScheme)
{
	ClearAllGradientDirections();
	std::stringstream ss(gradientScheme);

	//-----------------------------------------------
	// Read the nominal b-value
	//-----------------------------------------------
	int nominalB0=0;
	ss >> nominalB0;
	if ( nominalB0==0 )
		throw itk::ExceptionObject(__FILE__,__LINE__,"Invalid gradient description. It should starts with the nominal B0", "");
	SetNominalB0(nominalB0);

	bool stop=false;
	int nbB0=0;

	std::cout<<"- Gradient scheme (Nominal b-value="<<nominalB0<<"):"<<std::endl;

	//-----------------------------------------------
	// Read the full gradient scheme description
	//-----------------------------------------------
	while ( !stop )
	{
		char Command;
		ss >> Command;
		
		//-----------------------------------------------
		// Read a shell description
		//-----------------------------------------------
		if ( Command=='S' ) 
		{
			int nbDir, BVal;
			ss >> nbDir;
			ss >> BVal;
			if ( nbDir==0 || BVal==0 )
				throw itk::ExceptionObject(__FILE__,__LINE__,"Invalid shell gradient description. It should be 'S [nbDir] [bVal]' ", "");

			std::cout<<"  - One shell, "<<nbDir<<"dirs, b="<<BVal<<"."<<std::endl;
			AddOneShell(nbDir, BVal);
		}
		//-----------------------------------------------
		// Read a cube corner description
		//-----------------------------------------------
		else if ( Command=='C' )
		{
			int Nex, Norm;
			ss >> Nex;
			ss >> Norm;
			if ( Nex==0 || (Norm!=2 && Norm!=3) )
				throw itk::ExceptionObject(__FILE__,__LINE__,"Invalid cube gradient description. It should be 'C [Nex] [Norm=2 or 3]' ", "");

			std::cout<<"  - Cube corners of norm "<<Norm<<" (Nex="<<Nex<<")."<<std::endl;
			if ( Norm==2 )
				AddCubeCorners( Nex, 0) ;
			else
				AddCubeCorners( 0, Nex ) ;
		}
		//-----------------------------------------------
		// Read a number of B=0 description
		//-----------------------------------------------
		else if ( Command=='0' )
		{
			int Nex;
			ss >> Nex;
			if ( Nex==0  )
				throw itk::ExceptionObject(__FILE__,__LINE__,"Invalid number of B=0 description. It should be '0 [Nex]' ", "");
			
			nbB0 += Nex;
		}

		if ( ss.eof() ) stop=true;
	}

	//-----------------------------------------------
	// If no B=0 were set, add by default 1
	//-----------------------------------------------
	// tmp to improve
	// copy the diff vectors in a temporary variable because InsertElement( 0, vect3d ); doesn't work
	std::vector<GradientDirectionType> tmpVectors(m_DiffusionVectors->Size());
	for ( int i=0; i<m_DiffusionVectors->Size(); i++ )
		tmpVectors[i]=m_DiffusionVectors->GetElement(i);
	
	m_DiffusionVectors->Initialize();

	if ( nbB0==0 ) nbB0=1;
	std::cout<<"  - Number of baseline: "<<nbB0<<"."<<std::endl;
	GradientDirectionType vect3d;
	vect3d[0] = vect3d[1] = vect3d[2] = 0;
	for ( int i=0; i<nbB0; i++ )
		m_DiffusionVectors->InsertElement( m_DiffusionVectors->Size(), vect3d );

	for ( int i=0; i<tmpVectors.size(); i++ )
		m_DiffusionVectors->InsertElement( m_DiffusionVectors->Size(), tmpVectors[i] );


}

/**********************************************************************************************//**
 * \fn	void DWIGradients::SetNominalB0(double _B0)
 *
 * \brief	Sets the nominal b-value (for the gradient scaling). 
 *
 * \author	Benoit Scherrer
 * \date	August 2010
 *
 * \exception	itk::ExceptionObject	Thrown when exception. 
 *
 * \param	_B0	The nominal b-value 
*************************************************************************************************/
void DWIGradients::SetNominalB0(double _B0)
{
	if ( m_DiffusionVectors->Size() != 0 )
		throw itk::ExceptionObject(__FILE__,__LINE__,"Cannot change the nominal B0 after a call to AddOneShell", "DWIGradients::SetNominalB0");

	m_NominalB0 = _B0;
}

/**********************************************************************************************//**
 * \fn	void DWIGradients::AddSpiralShell(int nbDirections, int NbBValues, double MinB,
 * 		double MaxB)
 *
 * \brief	Adds a spiral shell. 
 *
 * \author	Benoit Scherrer
 * \date	August 2010
 *
 * \param	nbDirections	The total nb of directions . 
 * \param	NbBValues		The total nb of different b-values. 
 * \param	MinB			The minimum b. 
 * \param	MaxB			The maximum b. 
*************************************************************************************************/
void DWIGradients::AddSpiralShell(int nbDirections, int NbBValues, double MinB, double MaxB)
{
	const double *p = GetGradientList(nbDirections);

	GradientDirectionType vect3d;
	for ( int i=0 ; i<nbDirections ; i++ )
	{
		vect3d[0] = *p; vect3d[1] = *(p+1); vect3d[2] = *(p+2);
		p += 3;

		double b = MaxB;
		if ( NbBValues > 1 )
		{
			int numB = (int) (std::floor( 0.5+(double) (i * (NbBValues-1))/((double)(nbDirections-1))));
			b = MinB + ((double)numB)/((double)(NbBValues-1)) * (MaxB-MinB);
		}
		if ( m_Verbose ) cout << "B-Value: "<<b<<endl;

		// Set the correct magnitude for the _B0 value
		vect3d = vect3d*(sqrt(b/m_NominalB0));

		m_DiffusionVectors->InsertElement( m_DiffusionVectors->Size(), vect3d );
	}
}

/**********************************************************************************************//**
 * \fn	void DWIGradients::AddOneShell(int nbDirections, double _B0)
 *
 * \brief	Adds one shell 
 *
 * \author	Benoit Scherrer
 * \date	August 2010
 *
 * \exception	itk::ExceptionObject	Thrown when exception. 
 *
 * \param	nbDirections	The nb of directions. 
 * \param	_B0				The b-value for that shell. 
*************************************************************************************************/
void DWIGradients::AddOneShell(int nbDirections, double _B0)
{
	if ( m_NominalB0 == 0 )
		throw itk::ExceptionObject(__FILE__,__LINE__,"The nominal B0 should be set before calling AddOneShell", "DWIGradients::AddOneShell");

	const double *p = GetGradientList(nbDirections);

	GradientDirectionType vect3d;
	for ( int i=0 ; i<nbDirections ; i++ )
	{
		vect3d[0] = *p; vect3d[1] = *(p+1); vect3d[2] = *(p+2);
		p += 3;

		// Set the correct magnitude for the _B0 value
		vect3d = vect3d*(sqrt(_B0/m_NominalB0));
		m_DiffusionVectors->InsertElement( m_DiffusionVectors->Size(), vect3d );
	}
	
}

/**********************************************************************************************//**
 * \fn	void DWIGradients::AddNullGradients(int nb)
 *
 * \brief	Adds null gradients. 
 *
 * \author	Benoit Scherrer
 * \date	August 2010
 *
 * \param	nb	The nb of B=0 gradients to add. 
*************************************************************************************************/
void DWIGradients::AddNullGradients(int nb)
{
	GradientDirectionType vect3d;
	vect3d[0] = vect3d[1] = vect3d[2] = 0;

	for ( int i=0 ; i<nb ; i++ )
	{
		m_DiffusionVectors->InsertElement( m_DiffusionVectors->Size(), vect3d );
	}
}

/**********************************************************************************************//**
 * \fn	void DWIGradients::AddCubeCorners(int nbRepeatCubeN2, int nbRepeatCubeN3)
 *
 * \brief	Adds cube corners 
 *
 * \author	Benoit Scherrer
 * \date	August 2010
 *
 * \param	nbRepeatCubeN2	The nb of repetition of the 2-norm gradients to add
 * \param	nbRepeatCubeN3	The nb of repetition of the 3-norm gradients to add 
*************************************************************************************************/
void DWIGradients::AddCubeCorners(int nbRepeatCubeN2, int nbRepeatCubeN3)
{
	const int nbOffsetN2 = 6;
	const int offsetN2_X[]={ 1,  0,  1, -1,  0, -1 };
	const int offsetN2_Y[]={ 1,  1,  0,  1, -1,  0 };
	const int offsetN2_Z[]={ 0,  1,  1,  0,  1,  1 };

	const int nbOffsetN3 = 4;
	const int offsetN3_X[]={ 1, -1,  1, -1 };
	const int offsetN3_Y[]={ 1,  1, -1, -1 };
	const int offsetN3_Z[]={ 1,  1,  1,  1 };

	if ( nbRepeatCubeN2>0 )
	{
		GradientDirectionType vect3d;

		float sym=1;
		for ( int n=0; n<nbRepeatCubeN2 ; n++ )
		{
			for ( int i=0 ; i<nbOffsetN2 ; i++ )
			{
				vect3d[0] = offsetN2_X[i]*sym; 
				vect3d[1] = offsetN2_Y[i]*sym; 
				vect3d[2] = offsetN2_Z[i]*sym; 

				m_DiffusionVectors->InsertElement( m_DiffusionVectors->Size(), vect3d );
			}
			sym *= -1;
		}
	}

	if (  nbRepeatCubeN3>0 )
	{
		GradientDirectionType vect3d;

		float sym=1;
		for ( int n=0; n<nbRepeatCubeN3 ; n++ )
		{
			for ( int i=0 ; i<nbOffsetN3 ; i++ )
			{
				vect3d[0] = offsetN3_X[i]*sym; 
				vect3d[1] = offsetN3_Y[i]*sym; 
				vect3d[2] = offsetN3_Z[i]*sym; 

				m_DiffusionVectors->InsertElement( m_DiffusionVectors->Size(), vect3d );
			}
			sym *= -1;
		}
	}
}

void DWIGradients::AddRandomGradientsOnCubeEdges(int nbGradients)
{
	const int offsetN2_X[]={ 1,  0,  1, -1,  0, -1 };
	const int offsetN2_Y[]={ 1,  1,  0,  1, -1,  0 };
	const int offsetN2_Z[]={ 0,  1,  1,  0,  1,  1 };

	GradientDirectionType vect3d;

	srand((unsigned int)time(NULL));
	for ( int g=0; g<nbGradients; g++ )
	{
		int edge = rand()%6;
		vect3d[0] = offsetN2_X[edge]; 
		vect3d[1] = offsetN2_Y[edge]; 
		vect3d[2] = offsetN2_Z[edge]; 

		float r = 2*((float)rand())/((float)RAND_MAX)-1;	// [-1 1]
		if ( vect3d[0]==0 ) vect3d[0]=r;
		else if ( vect3d[1]==0 ) vect3d[1]=r;
		else vect3d[2]=r;

		m_DiffusionVectors->InsertElement( m_DiffusionVectors->Size(), vect3d );
	}
}

const double* DWIGradients::GetGradientList(int nbDirections)
{
	const double *p = NULL;
	switch ( nbDirections )
	{	
	case 3:
		p = Elec003;
		break;
	case 4:
		p = Elec004;
		break;
	case 5:
		p = Elec005;
		break;
	case 6:
		p = Elec006;
		break;
	case 7:
		p = Elec007;
		break;
	case 8:
		p = Elec008;
		break;
	case 9:
		p = Elec009;
		break;
	case 10:
		p = Elec010;
		break;
	case 11:
		p = Elec011;
		break;
	case 12:
		p = Elec012;
		break;
	case 13:
		p = Elec013;
		break;
	case 14:
		p = Elec014;
		break;
	case 15:
		p = Elec015;
		break;
	case 16:
		p = Elec016;
		break;
	case 17:
		p = Elec017;
		break;
	case 18:
		p = Elec018;
		break;
	case 19:
		p = Elec019;
		break;
	case 20:
		p = Elec020;
		break;
	case 21:
		p = Elec021;
		break;
	case 22:
		p = Elec022;
		break;
	case 23:
		p = Elec023;
		break;
	case 24:
		p = Elec024;
		break;
	case 25:
		p = Elec025;
		break;
	case 26:
		p = Elec026;
		break;
	case 27:
		p = Elec027;
		break;
	case 28:
		p = Elec028;
		break;
	case 29:
		p = Elec029;
		break;
	case 30:
		p = Elec030;
		break;
	case 31:
		p = Elec031;
		break;
	case 32:
		p = Elec032;
		break;
	case 33:
		p = Elec033;
		break;
	case 34:
		p = Elec034;
		break;
	case 35:
		p = Elec035;
		break;
	case 36:
		p = Elec036;
		break;
	case 37:
		p = Elec037;
		break;
	case 38:
		p = Elec038;
		break;
	case 39:
		p = Elec039;
		break;
	case 40:
		p = Elec040;
		break;
	case 41:
		p = Elec041;
		break;
	case 42:
		p = Elec042;
		break;
	case 43:
		p = Elec043;
		break;
	case 44:
		p = Elec044;
		break;
	case 45:
		p = Elec045;
		break;
	case 46:
		p = Elec046;
		break;
	case 47:
		p = Elec047;
		break;
	case 48:
		p = Elec048;
		break;
	case 49:
		p = Elec049;
		break;
	case 50:
		p = Elec050;
		break;
	case 51:
		p = Elec051;
		break;
	case 52:
		p = Elec052;
		break;
	case 53:
		p = Elec053;
		break;
	case 54:
		p = Elec054;
		break;
	case 55:
		p = Elec055;
		break;
	case 56:
		p = Elec056;
		break;
	case 57:
		p = Elec057;
		break;
	case 58:
		p = Elec058;
		break;
	case 59:
		p = Elec059;
		break;
	case 60:
		p = Elec060;
		break;
	case 61:
		p = Elec061;
		break;
	case 62:
		p = Elec062;
		break;
	case 63:
		p = Elec063;
		break;
	case 64:
		p = Elec064;
		break;
	case 65:
		p = Elec065;
		break;
	case 66:
		p = Elec066;
		break;
	case 67:
		p = Elec067;
		break;
	case 68:
		p = Elec068;
		break;
	case 69:
		p = Elec069;
		break;
	case 70:
		p = Elec070;
		break;
	case 71:
		p = Elec071;
		break;
	case 72:
		p = Elec072;
		break;
	case 73:
		p = Elec073;
		break;
	case 74:
		p = Elec074;
		break;
	case 75:
		p = Elec075;
		break;
	case 76:
		p = Elec076;
		break;
	case 77:
		p = Elec077;
		break;
	case 78:
		p = Elec078;
		break;
	case 79:
		p = Elec079;
		break;
	case 80:
		p = Elec080;
		break;
	case 81:
		p = Elec081;
		break;
	case 82:
		p = Elec082;
		break;
	case 83:
		p = Elec083;
		break;
	case 84:
		p = Elec084;
		break;
	case 85:
		p = Elec085;
		break;
	case 86:
		p = Elec086;
		break;
	case 87:
		p = Elec087;
		break;
	case 88:
		p = Elec088;
		break;
	case 89:
		p = Elec089;
		break;
	case 90:
		p = Elec090;
		break;
	case 91:
		p = Elec091;
		break;
	case 92:
		p = Elec092;
		break;
	case 93:
		p = Elec093;
		break;
	case 94:
		p = Elec094;
		break;
	case 95:
		p = Elec095;
		break;
	case 96:
		p = Elec096;
		break;
	case 97:
		p = Elec097;
		break;
	case 98:
		p = Elec098;
		break;
	case 99:
		p = Elec099;
		break;
	case 100:
		p = Elec100;
		break;
	case 101:
		p = Elec101;
		break;
	case 102:
		p = Elec102;
		break;
	case 103:
		p = Elec103;
		break;
	case 104:
		p = Elec104;
		break;
	case 105:
		p = Elec105;
		break;
	case 106:
		p = Elec106;
		break;
	case 107:
		p = Elec107;
		break;
	case 108:
		p = Elec108;
		break;
	case 109:
		p = Elec109;
		break;
	case 110:
		p = Elec110;
		break;
	case 111:
		p = Elec111;
		break;
	case 112:
		p = Elec112;
		break;
	case 113:
		p = Elec113;
		break;
	case 114:
		p = Elec114;
		break;
	case 115:
		p = Elec115;
		break;
	case 116:
		p = Elec116;
		break;
	case 117:
		p = Elec117;
		break;
	case 118:
		p = Elec118;
		break;
	case 119:
		p = Elec119;
		break;
	case 120:
		p = Elec120;
		break;
	case 121:
		p = Elec121;
		break;
	case 122:
		p = Elec122;
		break;
	case 123:
		p = Elec123;
		break;
	case 124:
		p = Elec124;
		break;
	case 125:
		p = Elec125;
		break;
	case 126:
		p = Elec126;
		break;
	case 127:
		p = Elec127;
		break;
	case 128:
		p = Elec128;
		break;
	case 129:
		p = Elec129;
		break;
	case 130:
		p = Elec130;
		break;
	case 131:
		p = Elec131;
		break;
	case 132:
		p = Elec132;
		break;
	case 133:
		p = Elec133;
		break;
	case 134:
		p = Elec134;
		break;
	case 135:
		p = Elec135;
		break;
	case 136:
		p = Elec136;
		break;
	case 137:
		p = Elec137;
		break;
	case 138:
		p = Elec138;
		break;
	case 139:
		p = Elec139;
		break;
	case 140:
		p = Elec140;
		break;
	case 141:
		p = Elec141;
		break;
	case 142:
		p = Elec142;
		break;
	case 143:
		p = Elec143;
		break;
	case 144:
		p = Elec144;
		break;
	case 145:
		p = Elec145;
		break;
	case 146:
		p = Elec146;
		break;
	case 147:
		p = Elec147;
		break;
	case 148:
		p = Elec148;
		break;
	case 149:
		p = Elec149;
		break;
	case 150:
		p = Elec150;
		break;
	case 151:
		p = Elec151;
		break;
	case 152:
		p = Elec152;
		break;
	case 153:
		p = Elec153;
		break;
	case 154:
		p = Elec154;
		break;
	case 155:
		p = Elec155;
		break;
	case 156:
		p = Elec156;
		break;
	case 157:
		p = Elec157;
		break;
	case 158:
		p = Elec158;
		break;
	case 159:
		p = Elec159;
		break;
	case 160:
		p = Elec160;
		break;
	case 161:
		p = Elec161;
		break;
	case 162:
		p = Elec162;
		break;
	case 163:
		p = Elec163;
		break;
	case 164:
		p = Elec164;
		break;
	case 165:
		p = Elec165;
		break;
	case 166:
		p = Elec166;
		break;
	case 167:
		p = Elec167;
		break;
	case 168:
		p = Elec168;
		break;
	case 169:
		p = Elec169;
		break;
	case 170:
		p = Elec170;
		break;
	case 171:
		p = Elec171;
		break;
	case 172:
		p = Elec172;
		break;
	case 173:
		p = Elec173;
		break;
	case 174:
		p = Elec174;
		break;
	case 175:
		p = Elec175;
		break;
	case 176:
		p = Elec176;
		break;
	case 177:
		p = Elec177;
		break;
	case 178:
		p = Elec178;
		break;
	case 179:
		p = Elec179;
		break;
	case 180:
		p = Elec180;
		break;
	case 181:
		p = Elec181;
		break;
	case 182:
		p = Elec182;
		break;
	case 183:
		p = Elec183;
		break;
	case 184:
		p = Elec184;
		break;
	case 185:
		p = Elec185;
		break;
	case 186:
		p = Elec186;
		break;
	case 187:
		p = Elec187;
		break;
	case 188:
		p = Elec188;
		break;
	case 189:
		p = Elec189;
		break;
	case 190:
		p = Elec190;
		break;
	case 191:
		p = Elec191;
		break;
	case 192:
		p = Elec192;
		break;
	case 193:
		p = Elec193;
		break;
	case 194:
		p = Elec194;
		break;
	case 195:
		p = Elec195;
		break;
	case 196:
		p = Elec196;
		break;
	case 197:
		p = Elec197;
		break;
	case 198:
		p = Elec198;
		break;
	case 199:
		p = Elec199;
		break;
	case 200:
		p = Elec200;
		break;
	case 201:
		p = Elec201;
		break;
	case 202:
		p = Elec202;
		break;
	case 203:
		p = Elec203;
		break;
	case 204:
		p = Elec204;
		break;
	case 205:
		p = Elec205;
		break;
	case 206:
		p = Elec206;
		break;
	case 207:
		p = Elec207;
		break;
	case 208:
		p = Elec208;
		break;
	case 209:
		p = Elec209;
		break;
	case 210:
		p = Elec210;
		break;
	case 211:
		p = Elec211;
		break;
	case 212:
		p = Elec212;
		break;
	case 213:
		p = Elec213;
		break;
	case 214:
		p = Elec214;
		break;
	case 215:
		p = Elec215;
		break;
	case 216:
		p = Elec216;
		break;
	case 217:
		p = Elec217;
		break;
	case 218:
		p = Elec218;
		break;
	case 219:
		p = Elec219;
		break;
	case 220:
		p = Elec220;
		break;
	case 221:
		p = Elec221;
		break;
	case 222:
		p = Elec222;
		break;
	case 223:
		p = Elec223;
		break;
	case 224:
		p = Elec224;
		break;
	case 225:
		p = Elec225;
		break;
	case 226:
		p = Elec226;
		break;
	case 227:
		p = Elec227;
		break;
	case 228:
		p = Elec228;
		break;
	case 229:
		p = Elec229;
		break;
	case 230:
		p = Elec230;
		break;
	case 231:
		p = Elec231;
		break;
	case 232:
		p = Elec232;
		break;
	case 233:
		p = Elec233;
		break;
	case 234:
		p = Elec234;
		break;
	case 235:
		p = Elec235;
		break;
	case 236:
		p = Elec236;
		break;
	case 237:
		p = Elec237;
		break;
	case 238:
		p = Elec238;
		break;
	case 239:
		p = Elec239;
		break;
	case 240:
		p = Elec240;
		break;
	case 241:
		p = Elec241;
		break;
	case 242:
		p = Elec242;
		break;
	case 243:
		p = Elec243;
		break;
	case 244:
		p = Elec244;
		break;
	case 245:
		p = Elec245;
		break;
	case 246:
		p = Elec246;
		break;
	case 247:
		p = Elec247;
		break;
	case 248:
		p = Elec248;
		break;
	case 249:
		p = Elec249;
		break;
	case 250:
		p = Elec250;
		break;
	case 251:
		p = Elec251;
		break;
	case 252:
		p = Elec252;
		break;
	case 253:
		p = Elec253;
		break;
	case 254:
		p = Elec254;
		break;
	case 255:
		p = Elec255;
		break;
	case 256:
		p = Elec256;
		break;
	case 755:
		p = Elec755;
		break;
	default:
		throw itk::ExceptionObject(__FILE__,__LINE__,"Invalid number of directions", "DWIGradients::AddOneShell");
	}


	return p;
}




} // end namespace DWI.
} // end namespace crl.

//***********************************************************************


