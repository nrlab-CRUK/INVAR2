/*
 * Functions to help pipeline development.
 */

/*
 * Log an exception to the logger as an error, including the stack trace.
 * Looks for InvocationTargetExceptions, which occur quite often, and logs
 * the cause of that exception, not the wrapper exception.
 */
def logException(e)
{
    def forLogging = e
    try
    {
        throw e
    }
    catch (java.lang.reflect.InvocationTargetException ite)
    {
        forLogging = e.targetException
    }
    catch (Throwable t)
    {
    }

    def sw = new StringWriter(1000)
    def pw = new PrintWriter(sw)
    forLogging.printStackTrace(pw)
    log.error sw.toString()
    throw e
}
