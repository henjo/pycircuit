#include "gnucap/u_lang.h"
#include "gnucap/c_comand.h"
#include "gnucap/globals.h"

#include <Python.h>

static char gnucap_doc[] = "Gnucap python interface module";

static PyObject *gnucap_command(PyObject *self, PyObject *args) {
  int ok;
  char *command;
  
  ok = PyArg_ParseTuple(args, "s", &command);
  
  if(ok) {
    SET_RUN_MODE xx(rBATCH);
    try {
      CMD::command(std::string(command), &CARD_LIST::card_list);
    } catch (Exception& e) {
      PyErr_SetString(PyExc_Exception, e.message().c_str());
      return NULL;
    }
  }
  
  Py_RETURN_NONE;
}

static PyMethodDef gnucap_methods[] = {
  {"command", gnucap_command, METH_VARARGS,
     "Run gnucap command"},
  {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC
initgnucap(void)
{
  Py_InitModule3("gnucap", gnucap_methods, 
                 "Gnucap python interface module");
}

