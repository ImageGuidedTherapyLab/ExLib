#include "itkCommand.h"

namespace crl
{

template < class TOptimizer >
class CommandIterationUpdate : public itk::Command
{
public:
  typedef  CommandIterationUpdate   Self;
  typedef  itk::Command             Superclass;
  typedef itk::SmartPointer<Self>  Pointer;
  itkNewMacro( Self );
protected:
  CommandIterationUpdate() {};
public:
  typedef TOptimizer     OptimizerType;
  typedef   const OptimizerType   *    OptimizerPointer;

  void Execute(itk::Object *caller, const itk::EventObject & event)
    {
      Execute( (const itk::Object *)caller, event);
    }

  void Execute(const itk::Object * object, const itk::EventObject & event)
    {
      OptimizerPointer optimizer =
        dynamic_cast< OptimizerPointer >( object );
      if( ! itk::IterationEvent().CheckEvent( &event ) )
        {
        return;
        }
      std::cout << optimizer->GetCurrentIteration() << "   ";
      std::cout << optimizer->GetValue();
      std::cout << std::endl;
    }
  
  
  // 
  // static function... you can do something like this:
  // crl::CommandIterationUpdate<OptimizerType>::CreateAndRegister(optimizer)
  // in your code.  You can get the observer as a return value if you like.
  //
  static Pointer CreateAndRegister( OptimizerType* optimizer )
    {
    Pointer observer = New();
    optimizer->AddObserver( itk::IterationEvent(), observer);
    return observer;
    }
  
};


};
