/**
 * Simple helper class for sending mail from Nextflow hooks.
 *
 * Requirements:
 *   - `mail` or `mailx` must be installed and configured on the system/container
 */
class Mailer {

    /**
     * Send a mail message.
     *
     * @param to Recipient address
     * @param subject Subject line
     * @param body Message body
     */
    static void send(String to, String subject, String body) {
        // Adjust command if your cluster uses `mailx` instead of `mail`
        def cmd = ['bash', '-lc', "mail -s \"${subject}\" ${to}"]
        def p = cmd.execute()

        // Write the body to stdin of the mail process
        p.outputStream.withWriter('UTF-8') { w -> w << body }
        p.outputStream.close()

        p.waitFor()

        if (p.exitValue() != 0) {
            def err = p.err.text.trim()
            println "[WARN] Mail send failed (exit ${p.exitValue()}): ${err}"
        }
    }
}
